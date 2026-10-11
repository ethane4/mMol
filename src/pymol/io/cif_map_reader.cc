/*
 * Read a molecule from CIF
 *
 * (c) 2014 Schrodinger, Inc.
 */

// Reading maps from mmCIF structure factors: symmetry expansion and Fourier
// transform.

#include "pymol/executive/executive.h"
#include "pocketfft_hdronly.h"
#include "pymol/io/cif_molecule_reader_internal.h"

constexpr auto EPSILON = std::numeric_limits<float>::epsilon();

/**
 * @brief Structure to hold reciprocal space parameters
 */
struct ReciprocalSpaceParams {
  glm::vec3 abc{};
  glm::vec3 angles{};
  float volume{};
};

/**
 * @brief Calculate the reciprocal space parameters for a given crystal
 * @param crystal The crystal object
 * @return The reciprocal space parameters
 * @note Fundamentals of Crystallography, 3rd edition (Table 2.1)
 */
ReciprocalSpaceParams CalculateReciprocalSpaceParam(const CCrystal& crystal) {
  auto dims = crystal.dims();
  auto a = dims[0];
  auto b = dims[1];
  auto c = dims[2];
  auto angles = crystal.angles();
  auto alpha = glm::radians(angles[0]);
  auto beta = glm::radians(angles[1]);
  auto gamma = glm::radians(angles[2]);
  auto cosAlpha = std::cos(alpha);
  auto cosBeta = std::cos(beta);
  auto cosGamma = std::cos(gamma);
  auto sinAlpha = std::sin(alpha);
  auto sinBeta = std::sin(beta);
  auto sinGamma = std::sin(gamma);
  auto volume =
      a * b * c *
      std::sqrt(1.0 - cosAlpha * cosAlpha - cosBeta * cosBeta -
                cosGamma * cosGamma + 2 * cosAlpha * cosBeta * cosGamma);
  auto inv_volume = 1.0 / volume;
  auto aStar = b * c * sinAlpha * inv_volume;
  auto bStar = a * c * sinBeta * inv_volume;
  auto cStar = a * b * sinGamma * inv_volume;
  auto cosAlphaStar = (cosBeta * cosGamma - cosAlpha) / (sinBeta * sinGamma);
  auto cosBetaStar = (cosAlpha * cosGamma - cosBeta) / (sinAlpha * sinGamma);
  auto cosGammaStar = (cosAlpha * cosBeta - cosGamma) / (sinAlpha * sinBeta);
  ReciprocalSpaceParams params{};
  params.abc = {aStar, bStar, cStar};
  params.angles = {cosAlphaStar, cosBetaStar, cosGammaStar};
  params.volume = static_cast<float>(volume);
  return params;
}

/**
 * @brief Calculate the squared distance in reciprocal space for a given hkl
 * @param param The reciprocal space parameters
 * @param hkl The Miller indices
 * @return The squared distance in reciprocal space
 * @note Fundamentals of Crystallography, 3rd edition (Eqn 2.17b)
 *      triclinic system is the most general case
 */
float CalculateReciprocalSpaceDistanceSquared(
    const ReciprocalSpaceParams& param, const glm::vec3& hkl) {
  auto a_star = param.abc[0];
  auto b_star = param.abc[1];
  auto c_star = param.abc[2];
  auto cos_alpha_star = param.angles[0];
  auto cos_beta_star = param.angles[1];
  auto cos_gamma_star = param.angles[2];
  auto h = hkl.x;
  auto k = hkl.y;
  auto l = hkl.z;
  return (h * h * a_star * a_star +                       //
          k * k * b_star * b_star +                       //
          l * l * c_star * c_star +                       //
          2 * h * k * a_star * b_star * cos_gamma_star +  //
          2 * h * l * a_star * c_star * cos_beta_star +   //
          2 * k * l * b_star * c_star * cos_alpha_star);
}

/**
 * @brief Structure to hold reflection information
 */
struct ReflectionInfo {
  glm::ivec3 hkl{};
  float fwt{};
  float phwt{};
  float fom{};
};

/**
 * @brief Class to handle the reflection data from a CIF data block
 */
class ReflnBlock {
 public:
  /**
   * @brief Constructor
   * @param datablock The CIF data block containing reflection data
   * @throws if any of the six required arrays are missing
   */
  ReflnBlock(const pymol::cif_data& datablock)
      : m_millerH(datablock.get_arr("_refln.index_h")),
        m_millerK(datablock.get_arr("_refln.index_k")),
        m_millerL(datablock.get_arr("_refln.index_l")),
        m_fwt(datablock.get_arr("_refln.pdbx_fwt")),
        m_phwt(datablock.get_arr("_refln.pdbx_phwt")),
        m_fom(datablock.get_arr("_refln.fom")) {
    if (!m_millerH || !m_millerK || !m_millerL || !m_fwt || !m_phwt || !m_fom) {
      m_size = 0;
      return;
    }
    m_size = m_millerH->size();
  }

  /**
   * @brief Get the reflection information for a given index
   * @param index The index of the reflection
   * @return The reflection information
   */
  pymol::Result<ReflectionInfo> getReflectionInfo(int index) const {
    if (index < 0 || index >= m_size) {
      return pymol::make_error("Index out of bounds");
    }

    return ReflectionInfo{
        glm::ivec3(m_millerH->as_i(index), m_millerK->as_i(index),
                   m_millerL->as_i(index)),
        m_fwt->as_f(index), m_phwt->as_f(index), m_fom->as_f(index)};
  }

  /**
   * @return The number of reflections in the block
   */
  auto size() const { return m_size; }

 private:
  const cif_array* m_millerH{};
  const cif_array* m_millerK{};
  const cif_array* m_millerL{};
  const cif_array* m_fwt{};
  const cif_array* m_phwt{};
  const cif_array* m_fom{};
  std::size_t m_size{};
};

/**
 * @brief Convert HKL to flat index
 * @param hkl The Miller indices
 * @param shape The shape of the map
 * @return The flat index corresponding to the HKL
 */
static std::size_t HKLToFlatIndex(const glm::ivec3& hkl,
                                  const glm::ivec3& shape) {
  return static_cast<std::size_t>(hkl.z + hkl.y * shape[2] +
                                  hkl.x * shape[1] * shape[2]);
}

/**
 * @brief Convert flat index to HKL
 * @param flat_idx The flat index
 * @param shape The shape of the map
 * @return The HKL corresponding to the flat index
 */
static glm::ivec3 FlatIndexToHKL(std::size_t flat_idx,
                                 const glm::ivec3& shape) {
  int iz = flat_idx % shape[2];
  int iy = (flat_idx / shape[2]) % shape[1];
  int ix = flat_idx / (shape[1] * shape[2]);
  return glm::ivec3(ix, iy, iz);
}

/**
 * @brief Create an ObjectMapState from a CFieldTyped map and symmetry
 * @param map The map to convert
 * @param sym The symmetry to use
 * @return The created ObjectMapState
 */
ObjectMapState ObjectMapStateFromField(PyMOLGlobals* G,
                                       const CFieldTyped<float>& map,
                                       std::unique_ptr<CSymmetry> sym) {
  ObjectMapState ms(G);
  glm::ivec3 grid_size(map.dim[0], map.dim[1], map.dim[2]);
  ms.MapSource = cMapSourceCCP4;
  ms.Active = true;
  ms.Symmetry.reset(sym.release());
  ms.FDim[0] = grid_size[0];
  ms.FDim[1] = grid_size[1];
  ms.FDim[2] = grid_size[2];
  ms.FDim[3] = 3;

  ms.Min[0] = 0;
  ms.Min[1] = 0;
  ms.Min[2] = 0;
  ms.Max[0] = grid_size.x - 1;
  ms.Max[1] = grid_size.y - 1;
  ms.Max[2] = grid_size.z - 1;

  ms.Div[0] = grid_size.x;
  ms.Div[1] = grid_size.y;
  ms.Div[2] = grid_size.z;

  ms.Field.reset(new Isofield(G, ms.FDim));

  float v[3];
  v[2] = (ms.Min[2]) / ((float)ms.Div[2]);
  v[1] = (ms.Min[1]) / ((float)ms.Div[1]);
  v[0] = (ms.Min[0]) / ((float)ms.Div[0]);

  transform33f3f(ms.Symmetry->Crystal.fracToReal(), v, ms.ExtentMin);

  v[2] = ((ms.FDim[2] - 1) + ms.Min[2]) / ((float)ms.Div[2]);
  v[1] = ((ms.FDim[1] - 1) + ms.Min[1]) / ((float)ms.Div[1]);
  v[0] = ((ms.FDim[0] - 1) + ms.Min[0]) / ((float)ms.Div[0]);

  transform33f3f(ms.Symmetry->Crystal.fracToReal(), v, ms.ExtentMax);

  auto& field = *ms.Field;
  assert(map.data.size() == field.data->data.size());
  std::copy(map.data.begin(), map.data.end(), field.data->data.begin());

  ms.Field->save_points = false;

  auto& crystal = ms.Symmetry->Crystal;

  ObjectMapStateRegeneratePoints(&ms);
  float frac_coords[3];
  float cart_coords[3];
  int corner_write_idx = 0;

  for (int z_loop = 0; z_loop < 2; ++z_loop) {
    int grid_idx_z = (z_loop == 0) ? ms.Min[2] : (ms.Min[2] + ms.FDim[2] - 1);
    frac_coords[2] =
        static_cast<float>(grid_idx_z) / static_cast<float>(ms.Div[2]);

    for (int y_loop = 0; y_loop < 2; ++y_loop) {
      int grid_idx_y = (y_loop == 0) ? ms.Min[1] : (ms.Min[1] + ms.FDim[1] - 1);
      frac_coords[1] =
          static_cast<float>(grid_idx_y) / static_cast<float>(ms.Div[1]);

      for (int x_loop = 0; x_loop < 2; ++x_loop) {
        int grid_idx_x =
            (x_loop == 0) ? ms.Min[0] : (ms.Min[0] + ms.FDim[0] - 1);
        frac_coords[0] =
            static_cast<float>(grid_idx_x) / static_cast<float>(ms.Div[0]);

        transform33f3f(ms.Symmetry->Crystal.fracToReal(), frac_coords,
                       cart_coords);

        ms.Corner[corner_write_idx * 3 + 0] = cart_coords[0];
        ms.Corner[corner_write_idx * 3 + 1] = cart_coords[1];
        ms.Corner[corner_write_idx * 3 + 2] = cart_coords[2];
        corner_write_idx++;
      }
    }
  }
  return ms;
}

/**
 * @brief Find the optimal grid size for the FFT
 * @param refln_block The reflection block containing reflection data
 * @param reciprocal_space_params The parameters for reciprocal space
 * @return The optimal grid size as a glm::ivec3
 */
glm::ivec3 FindOptimalGridSize(
    const ReflnBlock& refln_block,
    const ReciprocalSpaceParams& reciprocal_space_params) {
  std::array<int, 3> dimensions{};
  float inv_dist_sq{};
  for (std::size_t i = 0; i < refln_block.size(); ++i) {
    auto refl_info = refln_block.getReflectionInfo(i);
    if (!refl_info) {
      continue;
    }

    auto cand_inv_dist_sq = CalculateReciprocalSpaceDistanceSquared(
        reciprocal_space_params, refl_info->hkl);
    inv_dist_sq = std::max(inv_dist_sq, cand_inv_dist_sq);
  }
  double inv_dist = std::sqrt(inv_dist_sq);
  // Nyquist
  double nyquist_sampling_rate = 2.0;  // Perhaps make into a setting?
  glm::dvec3 cell_size{};
  for (int i = 0; i < 3; ++i) {
    cell_size[i] = std::max(cell_size[i], nyquist_sampling_rate * inv_dist /
                                              reciprocal_space_params.abc[i]);
  }

  glm::ivec3 grid_size{};
  for (int i = 0; i < 3; ++i) {
    auto round_up = static_cast<int>(std::ceil(cell_size[i]));
    grid_size[i] = pocketfft::detail::util::good_size_real(round_up);
  }
  return grid_size;
}

/**
 * @brief Calculate the phase shift for a given HKL and symmetry operation
 * @param hkl_in The input HKL indices
 * @param sym_op_T_3x1 The symmetry operation translation vector
 * @return The calculated phase shift (2 * PI * (h.trans))
 */
double CalculatePhaseShift(const glm::ivec3& hkl, const glm::vec3& trans) {
  double dot_product = glm::dot(glm::dvec3(hkl), glm::dvec3(trans));
  return 2.0 * glm::pi<double>() * dot_product;
}

/**
 * @brief Normalizes the complex map to float density map
 * @param reciprocal_space_params The parameters for reciprocal space
 * @param recip_grid The reciprocal grid data
 * @return Density Map as Float Field
 */
static CFieldTyped<float> NormalizeComplexMapToFloat(
    PyMOLGlobals* G, const ReciprocalSpaceParams& reciprocal_space_params,
    const CFieldTyped<std::complex<float>>& recip_grid) {
  CFieldTyped<float> map(reinterpret_cast<const int*>(recip_grid.dim.data()),
                         3);
  auto total_map_elements = map.size() / sizeof(float);

  float inv_Vcell = 1.0f;
  if (std::abs(reciprocal_space_params.volume) > EPSILON) {
    inv_Vcell = 1.0f / static_cast<float>(reciprocal_space_params.volume);
  }

  auto* flat_recip_data = recip_grid.ptr(0, 0, 0);
  auto* map_flat_data = map.ptr(0, 0, 0);

  if (!SettingGet<bool>(G, cSetting_normalize_ccp4_maps)) {
    for (std::size_t i = 0; i < total_map_elements; ++i) {
      map_flat_data[i] = flat_recip_data[i].real() * inv_Vcell;
    }
    return map;
  }

  double sum_rho_scaled = 0.0;
  double sum_sq_rho_scaled = 0.0;
  for (std::size_t i = 0; i < total_map_elements; ++i) {
    map_flat_data[i] = flat_recip_data[i].real() * inv_Vcell;
    sum_rho_scaled += map_flat_data[i];
    sum_sq_rho_scaled +=
        static_cast<double>(map_flat_data[i]) * map_flat_data[i];
  }

  // Normalize the map to mean ~0 and standard deviation ~1
  float current_mean = 0.0f;
  float current_sd = 0.0f;

  if (total_map_elements > 0) {
    current_mean = static_cast<float>(sum_rho_scaled / total_map_elements);
    double variance = (sum_sq_rho_scaled / total_map_elements) -
                      (static_cast<double>(current_mean) * current_mean);
    current_sd =
        (variance > 0) ? static_cast<float>(std::sqrt(variance)) : 0.0f;
  }

  if (std::abs(current_sd) > EPSILON) {  // Avoid division by zero if SD is tiny
    for (std::size_t i = 0; i < total_map_elements; ++i) {
      map_flat_data[i] = (map_flat_data[i] - current_mean) / current_sd;
    }
  } else {
    // Shift everything to force mean to zero
    if (std::abs(current_mean) > EPSILON) {
      for (std::size_t i = 0; i < total_map_elements; ++i) {
        map_flat_data[i] = map_flat_data[i] - current_mean;
      }
    }
  }
  return map;
}

/**
 * @brief Structure to hold symmetry matrices in GLM format
 */
struct SymMatsGLM {
  glm::mat3 rot;
  glm::vec3 trans;
};

/**
 * @brief Calculate the symmetry matrices in GLM format.
 * @param symmetry The symmetry object containing the symmetry matrices.
 * @return A vector of SymMatsGLM structures containing the rotation and
 * translation components of the symmetry matrices.
 */
static std::vector<SymMatsGLM> CalculateSymMatricesGLM(
    const CSymmetry& symmetry) {
  std::vector<SymMatsGLM> sym_matrices(symmetry.getNSymMat());
  for (int i = 0; i < symmetry.getNSymMat(); ++i) {
    auto* sym_matrix_4x4 = symmetry.getSymMat(i);
    const float R_j[9] = {
        sym_matrix_4x4[0], sym_matrix_4x4[1], sym_matrix_4x4[2],
        sym_matrix_4x4[4], sym_matrix_4x4[5], sym_matrix_4x4[6],
        sym_matrix_4x4[8], sym_matrix_4x4[9], sym_matrix_4x4[10]};
    sym_matrices[i].rot = glm::transpose(glm::make_mat3(R_j));
    sym_matrices[i].trans =
        glm::vec3(sym_matrix_4x4[3], sym_matrix_4x4[7], sym_matrix_4x4[11]);
  }
  return sym_matrices;
}

/**
 * @brief Checks if HKL indices are within the bounds of a grid.
 * @param idx The HKL indices to check.
 * @param shape The shape of the grid.
 * @return True if the indices are within bounds, false otherwise.
 */
static bool IsHKLWithinShape(const glm::ivec3& idx, const glm::ivec3& shape) {
  return (idx.x >= 0 && idx.x < shape[0] &&  //
          idx.y >= 0 && idx.y < shape[1] &&  //
          idx.z >= 0 && idx.z < shape[2]);
}

/**
 * @brief Apply symmetry matrices to the reflections and store the results in
 * the reciprocal space data.
 * @param sym_matrices The symmetry matrices to apply.
 * @param hkl_orig The original Miller indices.
 * @param shape The shape of the reciprocal space grid.
 * @param f_amp_weighted The weighted amplitude of the reflection.
 * @param phase_rad_orig The original phase in radians.
 * @param recip_1D_total_elements The total number of elements in the 1D
 * reciprocal space data.
 * @param[out] flat_recip_data The flat reciprocal space data array.
 * @param[out] is_grid_point_set A vector indicating whether a grid point has
 * been set.
 */
void ApplySymMatricesToReflections(
    const std::vector<SymMatsGLM>& sym_matrices, const glm::ivec3& hkl_orig,
    const glm::ivec3& shape, float f_amp_weighted, float phase_rad_orig,
    int recip_1D_total_elements, std::complex<float>* flat_recip_data,
    std::vector<std::uint8_t>& is_grid_point_set) {
  for (auto& sym_mat : sym_matrices) {
    glm::ivec3 hkl_symm = glm::round(sym_mat.rot * hkl_orig);
    auto phase_shift_rad = CalculatePhaseShift(hkl_orig, sym_mat.trans);

    auto F_symm_complex = std::polar(
        f_amp_weighted, phase_rad_orig + static_cast<float>(phase_shift_rad));

    int idx_h_s = (hkl_symm.x >= 0) ? hkl_symm.x : hkl_symm.x + shape[0];
    int idx_k_s = (hkl_symm.y >= 0) ? hkl_symm.y : hkl_symm.y + shape[1];
    int idx_l_s = (hkl_symm.z >= 0) ? hkl_symm.z : hkl_symm.z + shape[2];
    auto idx_hkl_s = glm::ivec3(idx_h_s, idx_k_s, idx_l_s);

    if (IsHKLWithinShape(idx_hkl_s, shape)) {
      auto flat_idx_s = HKLToFlatIndex(idx_hkl_s, shape);
      if (flat_idx_s < recip_1D_total_elements &&
          !is_grid_point_set[flat_idx_s]) {
        flat_recip_data[flat_idx_s] = F_symm_complex;
        is_grid_point_set[flat_idx_s] = true;
      }
    }
  }
}

/**
 * @brief Calculate inverse Friedel mate indices for a given HKL and shape.
 * @param hkl_orig The original HKL indices.
 * @param shape The shape of the reciprocal space grid.
 * @return The Friedel mate indices as a glm::ivec3.
 */
glm::ivec3 CalculateFriedelMateIndex(const glm::ivec3& hkl_orig,
                                     const glm::ivec3& shape) {
  auto h_idx = hkl_orig.x;
  auto k_idx = hkl_orig.y;
  auto l_idx = hkl_orig.z;

  int h = (h_idx < shape[0] / 2 + shape[0] % 2) ? h_idx : h_idx - shape[0];
  int k = (k_idx < shape[1] / 2 + shape[1] % 2) ? k_idx : k_idx - shape[1];
  int l = (l_idx < shape[2] / 2 + shape[2] % 2) ? l_idx : l_idx - shape[2];

  // Friedel mate indices
  int h_bar = -h;
  int k_bar = -k;
  int l_bar = -l;

  // Map (-h,-k,-l) to Friedel grid indices
  int idx_h_bar = (h_bar >= 0) ? h_bar : h_bar + shape[0];
  int idx_k_bar = (k_bar >= 0) ? k_bar : k_bar + shape[1];
  int idx_l_bar = (l_bar >= 0) ? l_bar : l_bar + shape[2];
  return glm::ivec3(idx_h_bar, idx_k_bar, idx_l_bar);
}

/**
 * @brief Add friedel pair to the reciprocal space data.
 * @param symmetry The symmetry object.
 * @param refln_block The reflection block containing reflection data.
 * @param reciprocal_space_params The parameters for reciprocal space.
 * @return A CFieldTyped<float> object representing the Fourier transformed map.
 */
void AddFriedelMate(const glm::ivec3& hkl_orig,
                    std::size_t recip_1D_total_elements,
                    const glm::ivec3& shape,
                    std::complex<float>* flat_recip_data,
                    std::vector<std::uint8_t>& is_grid_point_set) {
  if (!IsHKLWithinShape(hkl_orig, shape)) return;

  auto flat_idx_hkl = HKLToFlatIndex(hkl_orig, shape);
  if (flat_idx_hkl >= recip_1D_total_elements) return;

  auto idx_bar = CalculateFriedelMateIndex(hkl_orig, shape);

  if (IsHKLWithinShape(idx_bar, shape)) {
    auto flat_idx_bar = HKLToFlatIndex(idx_bar, shape);

    if (flat_idx_bar < recip_1D_total_elements &&
        !is_grid_point_set[flat_idx_bar]) {
      flat_recip_data[flat_idx_bar] = std::conj(flat_recip_data[flat_idx_hkl]);
      is_grid_point_set[flat_idx_bar] = true;
    }
    // If flat_idx_hkl == flat_idx_bar (centrosymmetric FFT point),
    // ensure it's real by zeroing negligible imaginary parts
    if (flat_idx_hkl == flat_idx_bar) {
      auto& recip_data = flat_recip_data[flat_idx_hkl];
      bool imag_is_negligible =
          std::abs(recip_data.imag()) < 1e-5f * std::abs(recip_data.real()) ||
          std::abs(recip_data.imag()) < 1e-5f;
      if (imag_is_negligible) {
        recip_data.imag(0.0f);
      }
    }
  }
}

/**
 * @brief Perform a Fast Fourier Transform in place on the reciprocal grid.
 * @param grid The reciprocal grid data.
 */
static void PerformFastFourierTransformInPlace(
    CFieldTyped<std::complex<float>>& grid) {
  auto* flat_grid = grid.ptr(0, 0, 0);
  pocketfft::shape_t pocket_shape{static_cast<std::size_t>(grid.dim[0]),
                                  static_cast<std::size_t>(grid.dim[1]),
                                  static_cast<std::size_t>(grid.dim[2])};
  auto base_size = grid.base_size;

  // CField is internally a 1D array, so we need to calculate the strides
  // for the 3D data. Z is contiguous; thus the axis order is ZYX.
  glm::uvec3 u_strides = {base_size * pocket_shape[2] * pocket_shape[1],
                          base_size * pocket_shape[2], base_size};
  pocketfft::stride_t stride_in(3);

  for (int i = 0; i < 3; ++i) {
    stride_in[i] = static_cast<pocketfft::stride_t::value_type>(u_strides[i]);
  }
  auto stride_out = stride_in;
  pocketfft::shape_t axes{2, 1, 0};  // ZYX
  bool direction = pocketfft::FORWARD;
  float norm = 1.0f;

  pocketfft::c2c<float>(pocket_shape, stride_in, stride_out, axes, direction,
                        flat_grid, flat_grid, norm);
}

/**
 * @brief Fourier transform structure factors to a map.
 * @param symmetry The symmetry object.
 * @param refln_block The reflection block containing reflection data.
 * @param reciprocal_space_params The parameters for reciprocal space.
 * @return A CFieldTyped<float> object representing the Fourier transformed map.
 */
CFieldTyped<float> FourierTransformStructureFactorsToMap(
    PyMOLGlobals* G, const CSymmetry& symmetry, const ReflnBlock& refln_block,
    const ReciprocalSpaceParams& reciprocal_space_params) {
  auto shape = FindOptimalGridSize(refln_block, reciprocal_space_params);
  CFieldTyped<std::complex<float>> recip_grid(glm::value_ptr(shape), 3);

  auto recip_1D_total_elements = shape[0] * shape[1] * shape[2];
  std::vector<std::uint8_t> is_grid_point_set(recip_1D_total_elements, 0);

  auto* flat_recip_data = recip_grid.ptr(0, 0, 0);

  auto sym_matrices = CalculateSymMatricesGLM(symmetry);

  for (std::size_t i = 0; i < refln_block.size(); ++i) {
    auto refl_info = refln_block.getReflectionInfo(i);
    if (!refl_info) {
      continue;
    }

    auto hkl_orig = refl_info->hkl;
    auto f_amp_raw = refl_info->fwt;
    auto ph_deg_orig = refl_info->phwt;
    auto fom_val = refl_info->fom;

    auto f_amp_weighted = f_amp_raw * fom_val;
    if (f_amp_weighted == 0.0f &&
        !(hkl_orig.x == 0 && hkl_orig.y == 0 &&
          hkl_orig.z == 0)) {  // Skip zero amplitude reflections unless F000
      continue;
    }

    auto phase_rad_orig = glm::radians(ph_deg_orig);

    ApplySymMatricesToReflections(sym_matrices, hkl_orig, shape, f_amp_weighted,
                                  phase_rad_orig, recip_1D_total_elements,
                                  flat_recip_data, is_grid_point_set);
  }  // End loop over reflections from CIF

  for (std::size_t h_idx = 0; h_idx < shape[0]; ++h_idx) {
    for (std::size_t k_idx = 0; k_idx < shape[1]; ++k_idx) {
      for (std::size_t l_idx = 0; l_idx < shape[2]; ++l_idx) {
        auto flat_idx_hkl =
            HKLToFlatIndex(glm::ivec3(h_idx, k_idx, l_idx), shape);
        if (is_grid_point_set[flat_idx_hkl]) {
          AddFriedelMate(glm::ivec3(h_idx, k_idx, l_idx),
                         recip_1D_total_elements, shape, flat_recip_data,
                         is_grid_point_set);
        }
      }
    }
  }

  PerformFastFourierTransformInPlace(recip_grid);
  auto map = NormalizeComplexMapToFloat(G, reciprocal_space_params, recip_grid);
  return map;
}

pymol::Result<ObjectMap*> ObjectMapReadCifStr(
    PyMOLGlobals* G, const pymol::cif_data& datablock) {
  // Validate required arrays before proceeding
  if (!datablock.get_arr("_refln.pdbx_fwt") ||
      !datablock.get_arr("_refln.pdbx_phwt")) {
    return pymol::make_error(
        "CIF reflection data missing required map coefficients "
        "(_refln.pdbx_FWT / _refln.pdbx_PHWT)");
  }

  auto I = std::make_unique<ObjectMap>(G);
  initializeTTT44f(I->TTT);
  I->TTTFlag = false;

  ReflnBlock refln_block(datablock);
  std::unique_ptr<CSymmetry> sym(read_symmetry(G, &datablock));
  if (!sym) {
    return pymol::make_error(
        "CIF reflection data missing symmetry/cell parameters");
  }

  auto reciprocal_space_params = CalculateReciprocalSpaceParam(sym->Crystal);
  auto map = FourierTransformStructureFactorsToMap(G, *sym, refln_block,
                                                   reciprocal_space_params);

  I->State.push_back(ObjectMapStateFromField(G, map, std::move(sym)));
  ObjectMapUpdateExtents(I.get());

  return I.release();
}
