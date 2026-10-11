/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2004 by Warren Lyford Delano of DeLano Scientific.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/

// Maps from selections: VDW masks, Gaussian densities and Coulomb potentials.

#include "pymol/base/feedback.h"
#include "pymol/select/selector.h"
#include "pymol/executive/executive.h"
#include "pymol/select/selector_internal.h"
int SelectorMapMaskVDW(PyMOLGlobals* G, int sele1, ObjectMapState* oMap,
                       float buffer, int state) {
  CSelector* I = G->Selector;
  float* v2;
  int n1;
  int a, b, c;
  int at;
  int s;
  ObjectMolecule* obj;
  CoordSet* cs;
  int state1, state2;
  int once_flag;

  c = 0;
  n1 = 0;
  SelectorUpdateTable(G, state, -1);

  const size_t table_size = I->Table.size();
  auto coords_flat = std::vector<float>(table_size * 3);
  auto* coords = pymol::reshape<3>(coords_flat.data());
  auto Flag1 = std::vector<MapFlag_t>(table_size, 0);

  for (a = 0; a < I->Table.size(); a++) {
    at = I->Table[a].atom;
    obj = I->Obj[I->Table[a].model];
    s = obj->AtomInfo[at].selEntry;
    if (SelectorIsMember(G, s, sele1)) {
      once_flag = true;
      for (state2 = 0; state2 < obj->NCSet; state2++) {
        if (state < 0) once_flag = false;
        if (!once_flag)
          state1 = state2;
        else
          state1 = state;
        if (state1 < obj->NCSet)
          cs = obj->CSet[state1];
        else
          cs = nullptr;
        if (cs) {
          if (CoordSetGetAtomVertex(cs, at, coords[a])) {
            Flag1[a] = true;
            n1++;
          }
        }
        if (once_flag) break;
      }
    }
  }
  /* now create and apply voxel map */
  c = 0;
  if (n1) {
    std::unique_ptr<MapType> map(new MapType(G, -(buffer + MAX_VDW),
                                             pymol::flatten(coords), table_size,
                                             nullptr, Flag1.data()));
    if (map) {
      for (a = oMap->Min[0]; a <= oMap->Max[0]; a++) {
        for (b = oMap->Min[1]; b <= oMap->Max[1]; b++) {
          for (c = oMap->Min[2]; c <= oMap->Max[2]; c++) {
            F3(oMap->Field->data, a, b, c) = 0.0;

            v2 = F4Ptr(oMap->Field->points, a, b, c, 0);

            for (const auto j : MapEIter(*map, v2)) {
              const auto* ai =
                  I->Obj[I->Table[j].model]->AtomInfo + I->Table[j].atom;
              if (within3f(coords[j], v2, ai->vdw + buffer)) {
                F3(oMap->Field->data, a, b, c) = 1.0;
              }
            }
          }
        }
      }
      oMap->Active = true;
    }
  }
  return (c);
}

static double max2d(double a, double b) {
  if (a > b)
    return a;
  else
    return b;
}

static double max6d(double a, double b, double c, double d, double e,
                    double f) {
  if (d > a) a = d;
  if (e > b) b = e;
  if (f > c) c = f;
  if (b > a) a = b;
  if (c > a) a = c;
  return (a);
}

#define D_SMALL10 1e-10

typedef double AtomSF[11];
int SelectorMapGaussian(PyMOLGlobals* G, int sele1, ObjectMapState* oMap,
                        float buffer, int state, int normalize, int use_max,
                        int quiet, float resolution) {
  CSelector* I = G->Selector;
  float* v2;
  int n1, n2;
  int a, b, c;
  int at;
  int s, idx;
  AtomInfoType* ai;
  ObjectMolecule* obj;
  CoordSet* cs;
  int state1, state2;
  float *point = nullptr, *fp;
  int *sfidx = nullptr, *ip;
  float *b_factor = nullptr, *bf, bfact;
  float *occup = nullptr, *oc;
  int prot;
  int once_flag;
  float d, e_val;
  double sum, sumsq;
  float mean, stdev;
  double sf[256][11], *sfp;
  AtomSF* atom_sf = nullptr;
  double b_adjust = (double)SettingGetGlobal_f(G, cSetting_gaussian_b_adjust);
  double elim = 7.0;
  double rcut2;
  float rcut;
  float max_rcut = 0.0F;
  float b_floor = SettingGetGlobal_f(G, cSetting_gaussian_b_floor);
  float blur_factor = 1.0F;

  {
    if (resolution < R_SMALL4)
      resolution = SettingGetGlobal_f(G, cSetting_gaussian_resolution);
    if (resolution < 1.0) resolution = 1.0F;
    blur_factor = 2.0F / resolution;
    /* a gaussian_resolution of 2.0 is considered perfect ? Hmm...where is this
     * from??? */
  }

  if (b_adjust > 500.0)
    b_adjust = 500.0; /* constrain to be somewhat reasonable */

  for (a = 0; a < 256; a++) {
    sf[a][0] = -1.0;
  }

  sf[cAN_H][0] = 0.493002;
  sf[cAN_H][1] = 10.510900;
  sf[cAN_H][2] = 0.322912;
  sf[cAN_H][3] = 26.125700;
  sf[cAN_H][4] = 0.140191;
  sf[cAN_H][5] = 3.142360;
  sf[cAN_H][6] = 0.040810;
  sf[cAN_H][7] = 57.799698;
  sf[cAN_H][8] = 0.003038;
  sf[cAN_H][9] = 0.0;

  /* LP currently using scattering factors of carbon
     (Roche Pocket viewer relies upon this behavior) */

  sf[cAN_LP][0] = 2.310000;
  sf[cAN_LP][1] = 20.843899;
  sf[cAN_LP][2] = 1.020000;
  sf[cAN_LP][3] = 10.207500;
  sf[cAN_LP][4] = 1.588600;
  sf[cAN_LP][5] = 0.568700;
  sf[cAN_LP][6] = 0.865000;
  sf[cAN_LP][7] = 51.651199;
  sf[cAN_LP][8] = 0.215600;
  sf[cAN_LP][9] = 0.0;

  sf[cAN_C][0] = 2.310000;
  sf[cAN_C][1] = 20.843899;
  sf[cAN_C][2] = 1.020000;
  sf[cAN_C][3] = 10.207500;
  sf[cAN_C][4] = 1.588600;
  sf[cAN_C][5] = 0.568700;
  sf[cAN_C][6] = 0.865000;
  sf[cAN_C][7] = 51.651199;
  sf[cAN_C][8] = 0.215600;
  sf[cAN_C][9] = 0.0;

  sf[cAN_O][0] = 3.048500;
  sf[cAN_O][1] = 13.277100;
  sf[cAN_O][2] = 2.286800;
  sf[cAN_O][3] = 5.701100;
  sf[cAN_O][4] = 1.546300;
  sf[cAN_O][5] = 0.323900;
  sf[cAN_O][6] = 0.867000;
  sf[cAN_O][7] = 32.908897;
  sf[cAN_O][8] = 0.250800;
  sf[cAN_O][9] = 0.0;

  sf[cAN_N][0] = 12.212600;
  sf[cAN_N][1] = 0.005700;
  sf[cAN_N][2] = 3.132200;
  sf[cAN_N][3] = 9.893300;
  sf[cAN_N][4] = 2.012500;
  sf[cAN_N][5] = 28.997499;
  sf[cAN_N][6] = 1.166300;
  sf[cAN_N][7] = 0.582600;
  sf[cAN_N][8] = -11.528999;
  sf[cAN_N][9] = 0.0;

  sf[cAN_S][0] = 6.905300;
  sf[cAN_S][1] = 1.467900;
  sf[cAN_S][2] = 5.203400;
  sf[cAN_S][3] = 22.215099;
  sf[cAN_S][4] = 1.437900;
  sf[cAN_S][5] = 0.253600;
  sf[cAN_S][6] = 1.586300;
  sf[cAN_S][7] = 56.172001;
  sf[cAN_S][8] = 0.866900;
  sf[cAN_S][9] = 0.0;

  sf[cAN_Cl][0] = 11.460400;
  sf[cAN_Cl][1] = 0.010400;
  sf[cAN_Cl][2] = 7.196400;
  sf[cAN_Cl][3] = 1.166200;
  sf[cAN_Cl][4] = 6.255600;
  sf[cAN_Cl][5] = 18.519400;
  sf[cAN_Cl][6] = 1.645500;
  sf[cAN_Cl][7] = 47.778400;
  sf[cAN_Cl][8] = 0.866900;
  sf[cAN_Cl][9] = 0.0;

  sf[cAN_Br][0] = 17.178900;
  sf[cAN_Br][1] = 2.172300;
  sf[cAN_Br][2] = 5.235800;
  sf[cAN_Br][3] = 16.579599;
  sf[cAN_Br][4] = 5.637700;
  sf[cAN_Br][5] = 0.260900;
  sf[cAN_Br][6] = 3.985100;
  sf[cAN_Br][7] = 41.432800;
  sf[cAN_Br][8] = 2.955700;
  sf[cAN_Br][9] = 0.0;

  sf[cAN_I][0] = 20.147200;
  sf[cAN_I][1] = 4.347000;
  sf[cAN_I][2] = 18.994900;
  sf[cAN_I][3] = 0.381400;
  sf[cAN_I][4] = 7.513800;
  sf[cAN_I][5] = 27.765999;
  sf[cAN_I][6] = 2.273500;
  sf[cAN_I][7] = 66.877602;
  sf[cAN_I][8] = 4.071200;
  sf[cAN_I][9] = 0.0;

  sf[cAN_F][0] = 3.539200;
  sf[cAN_F][1] = 10.282499;
  sf[cAN_F][2] = 2.641200;
  sf[cAN_F][3] = 4.294400;
  sf[cAN_F][4] = 1.517000;
  sf[cAN_F][5] = 0.261500;
  sf[cAN_F][6] = 1.024300;
  sf[cAN_F][7] = 26.147600;
  sf[cAN_F][8] = 0.277600;
  sf[cAN_F][9] = 0.0;

  sf[cAN_K][0] = 8.218599;
  sf[cAN_K][1] = 12.794900;
  sf[cAN_K][2] = 7.439800;
  sf[cAN_K][3] = 0.774800;
  sf[cAN_K][4] = 1.051900;
  sf[cAN_K][5] = 213.186996;
  sf[cAN_K][6] = 0.865900;
  sf[cAN_K][7] = 41.684097;
  sf[cAN_K][8] = 1.422800;
  sf[cAN_K][9] = 0.0;

  sf[cAN_Mg][0] = 5.420400;
  sf[cAN_Mg][1] = 2.827500;
  sf[cAN_Mg][2] = 2.173500;
  sf[cAN_Mg][3] = 79.261101;
  sf[cAN_Mg][4] = 1.226900;
  sf[cAN_Mg][5] = 0.380800;
  sf[cAN_Mg][6] = 2.307300;
  sf[cAN_Mg][7] = 7.193700;
  sf[cAN_Mg][8] = 0.858400;
  sf[cAN_Mg][9] = 0.0;

  sf[cAN_Na][0] = 4.762600;
  sf[cAN_Na][1] = 3.285000;
  sf[cAN_Na][2] = 3.173600;
  sf[cAN_Na][3] = 8.842199;
  sf[cAN_Na][4] = 1.267400;
  sf[cAN_Na][5] = 0.313600;
  sf[cAN_Na][6] = 1.112800;
  sf[cAN_Na][7] = 129.423996;
  sf[cAN_Na][8] = 0.676000;
  sf[cAN_Na][9] = 0.0;

  sf[cAN_P][0] = 6.434500;
  sf[cAN_P][1] = 1.906700;
  sf[cAN_P][2] = 4.179100;
  sf[cAN_P][3] = 27.157000;
  sf[cAN_P][4] = 1.780000;
  sf[cAN_P][5] = 0.526000;
  sf[cAN_P][6] = 1.490800;
  sf[cAN_P][7] = 68.164497;
  sf[cAN_P][8] = 1.114900;
  sf[cAN_P][9] = 0.0;

  sf[cAN_Zn][0] = 14.074300;
  sf[cAN_Zn][1] = 3.265500;
  sf[cAN_Zn][2] = 7.031800;
  sf[cAN_Zn][3] = 0.233300;
  sf[cAN_Zn][4] = 5.162500;
  sf[cAN_Zn][5] = 10.316299;
  sf[cAN_Zn][6] = 2.410000;
  sf[cAN_Zn][7] = 58.709702;
  sf[cAN_Zn][8] = 1.304100;
  sf[cAN_Zn][9] = 0.0;

  sf[cAN_Ca][0] = 8.626600;
  sf[cAN_Ca][1] = 10.442100;
  sf[cAN_Ca][2] = 7.387300;
  sf[cAN_Ca][3] = 0.659900;
  sf[cAN_Ca][4] = 1.589900;
  sf[cAN_Ca][5] = 85.748398;
  sf[cAN_Ca][6] = 1.021100;
  sf[cAN_Ca][7] = 178.436996;
  sf[cAN_Ca][8] = 1.375100;
  sf[cAN_Ca][9] = 0.0;

  sf[cAN_Cu][0] = 13.337999;
  sf[cAN_Cu][1] = 3.582800;
  sf[cAN_Cu][2] = 7.167600;
  sf[cAN_Cu][3] = 0.247000;
  sf[cAN_Cu][4] = 5.615800;
  sf[cAN_Cu][5] = 11.396600;
  sf[cAN_Cu][6] = 1.673500;
  sf[cAN_Cu][7] = 64.812599;
  sf[cAN_Cu][8] = 1.191000;
  sf[cAN_Cu][9] = 0.0;

  sf[cAN_Fe][0] = 11.769500;
  sf[cAN_Fe][1] = 4.761100;
  sf[cAN_Fe][2] = 7.357300;
  sf[cAN_Fe][3] = 0.307200;
  sf[cAN_Fe][4] = 3.522200;
  sf[cAN_Fe][5] = 15.353500;
  sf[cAN_Fe][6] = 2.304500;
  sf[cAN_Fe][7] = 76.880501;
  sf[cAN_Fe][8] = 1.036900;
  sf[cAN_Fe][9] = 0.0;

  sf[cAN_Se][0] = 17.000599;
  sf[cAN_Se][1] = 2.409800;
  sf[cAN_Se][2] = 5.819600;
  sf[cAN_Se][3] = 0.272600;
  sf[cAN_Se][4] = 3.973100;
  sf[cAN_Se][5] = 15.237200;
  sf[cAN_Se][6] = 4.354300;
  sf[cAN_Se][7] = 43.816299;
  sf[cAN_Se][8] = 2.840900;
  sf[cAN_Se][9] = 0.0;

  buffer += MAX_VDW;
  c = 0;
  n1 = 0;
  if (state >= cSelectorUpdateTableEffectiveStates) {
    SelectorUpdateTable(G, state, -1);
  } else {
    PRINTFB(G, FB_ObjectMap, FB_Warnings)
    " %s-Warning: state = %d\n", __func__, state ENDFB(G);
    assert(false);  // no mercy for debug build
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  }
  for (a = 0; a < I->Table.size(); a++) {
    at = I->Table[a].atom;
    obj = I->Obj[I->Table[a].model];
    s = obj->AtomInfo[at].selEntry;
    if (SelectorIsMember(G, s, sele1)) {
      once_flag = true;
      for (state1 = 0; state1 < obj->NCSet; state1++) {
        if (state < 0) once_flag = false;
        if (!once_flag)
          state2 = state1;
        else
          state2 = state;
        if (state2 < obj->NCSet)
          cs = obj->CSet[state2];
        else
          cs = nullptr;
        if (cs) {
          idx = cs->atmToIdx(at);
          if (idx >= 0) {
            n1++;
          }
        }
        if (once_flag) break;
      }
    }
  }
  point = pymol::malloc<float>(3 * n1);
  sfidx = pymol::malloc<int>(n1);
  b_factor = pymol::malloc<float>(n1);
  occup = pymol::malloc<float>(n1);
  atom_sf = pymol::malloc<AtomSF>(n1);

  if (!quiet) {
    PRINTFB(G, FB_ObjectMap, FB_Details)
    " ObjectMap: Computing Gaussian map for %d atom positions.\n", n1 ENDFB(G);
  }

  n1 = 0;
  fp = point;
  ip = sfidx;
  bf = b_factor;
  oc = occup;
  for (a = 0; a < I->Table.size(); a++) {
    at = I->Table[a].atom;
    obj = I->Obj[I->Table[a].model];
    ai = obj->AtomInfo + at;
    s = ai->selEntry;
    if (SelectorIsMember(G, s, sele1)) {
      once_flag = true;
      for (state1 = 0; state1 < obj->NCSet; state1++) {
        if (state < 0) once_flag = false;
        if (!once_flag)
          state2 = state1;
        else
          state2 = state;
        if (state2 < obj->NCSet)
          cs = obj->CSet[state2];
        else
          cs = nullptr;
        if (cs) {
          if (CoordSetGetAtomVertex(cs, at, fp)) {
            prot = ai->protons;
            if (sf[prot][0] == -1.0F) prot = cAN_C;
            bfact = ai->b + (float)b_adjust;
            if (bfact < b_floor) bfact = b_floor;
            if ((bfact > R_SMALL4) && (ai->q > R_SMALL4)) {
              fp += 3;
              *(ip++) = prot;
              *(bf++) = bfact;
              *(oc++) = ai->q;
              n1++;
            }
          }
        }
        if (once_flag) break;
      }
    }
  }

  for (a = 0; a < n1; a++) {
    double* src_sf;

    src_sf = &sf[sfidx[a]][0];
    bfact = b_factor[a];

    for (b = 0; b < 10; b += 2) {
      double sfa, sfb;
      sfa = src_sf[b];
      sfb = src_sf[b + 1];

      atom_sf[a][b] = occup[a] * sfa * pow(sqrt1d(4 * PI / (sfb + bfact)), 3.0);
      atom_sf[a][b + 1] = 4 * PI * PI / (sfb + bfact);
    }

    rcut2 = max6d(
        0.0,
        (elim + log(max2d(fabs(atom_sf[a][0]), D_SMALL10))) / atom_sf[a][1],
        (elim + log(max2d(fabs(atom_sf[a][2]), D_SMALL10))) / atom_sf[a][3],
        (elim + log(max2d(fabs(atom_sf[a][4]), D_SMALL10))) / atom_sf[a][5],
        (elim + log(max2d(fabs(atom_sf[a][6]), D_SMALL10))) / atom_sf[a][7],
        (elim + log(max2d(fabs(atom_sf[a][8]), D_SMALL10))) / atom_sf[a][9]);
    rcut = ((float)sqrt1d(rcut2)) / blur_factor;
    atom_sf[a][10] = rcut;
    if (max_rcut < rcut) max_rcut = rcut;
  }

  /* now create and apply voxel map */
  c = 0;
  if (n1) {
    n2 = 0;
    std::unique_ptr<MapType> map(new MapType(G, -max_rcut, point, n1, nullptr));
    if (map) {
      sum = 0.0;
      sumsq = 0.0;
      for (a = oMap->Min[0]; a <= oMap->Max[0]; a++) {
        OrthoBusyFast(G, a - oMap->Min[0], oMap->Max[0] - oMap->Min[0] + 1);
        for (b = oMap->Min[1]; b <= oMap->Max[1]; b++) {
          for (c = oMap->Min[2]; c <= oMap->Max[2]; c++) {
            e_val = 0.0;
            v2 = F4Ptr(oMap->Field->points, a, b, c, 0);
            if (use_max) {
              float e_partial;
              for (const auto j : MapEIter(*map, v2)) {
                d = (float)diff3f(point + 3 * j, v2) * blur_factor;
                /* scale up width */
                sfp = atom_sf[j];
                if (d < sfp[10]) {
                  d = d * d;
                  if (d < R_SMALL8) d = R_SMALL8;
                  e_partial = (float)((sfp[0] * exp(-sfp[1] * d)) +
                                      (sfp[2] * exp(-sfp[3] * d)) +
                                      (sfp[4] * exp(-sfp[5] * d)) +
                                      (sfp[6] * exp(-sfp[7] * d)) +
                                      (sfp[8] * exp(-sfp[9] * d))) *
                              blur_factor;
                  /* scale down intensity */
                  if (e_partial > e_val) e_val = e_partial;
                }
              }
            } else {
              for (const auto j : MapEIter(*map, v2)) {
                d = (float)diff3f(point + 3 * j, v2) * blur_factor;
                /* scale up width */
                sfp = atom_sf[j];
                if (d < sfp[10]) {
                  d = d * d;
                  if (d < R_SMALL8) d = R_SMALL8;
                  e_val += (float)((sfp[0] * exp(-sfp[1] * d)) +
                                   (sfp[2] * exp(-sfp[3] * d)) +
                                   (sfp[4] * exp(-sfp[5] * d)) +
                                   (sfp[6] * exp(-sfp[7] * d)) +
                                   (sfp[8] * exp(-sfp[9] * d))) *
                           blur_factor;
                  /* scale down intensity */
                }
              }
            }
            F3(oMap->Field->data, a, b, c) = e_val;
            sum += e_val;
            sumsq += (e_val * e_val);
            n2++;
          }
        }
      }
      mean = (float)(sum / n2);
      stdev = (float)sqrt1d((sumsq - (sum * sum / n2)) / (n2 - 1));
      if (normalize) {
        if (!quiet) {
          PRINTFB(G, FB_ObjectMap, FB_Details)
          " ObjectMap: Normalizing: mean = %8.6f & stdev = %8.6f.\n", mean,
              stdev ENDFB(G);
        }

        if (stdev < R_SMALL8) stdev = R_SMALL8;

        for (a = oMap->Min[0]; a <= oMap->Max[0]; a++) {
          for (b = oMap->Min[1]; b <= oMap->Max[1]; b++) {
            for (c = oMap->Min[2]; c <= oMap->Max[2]; c++) {
              fp = F3Ptr(oMap->Field->data, a, b, c);

              *fp = (*fp - mean) / stdev;
            }
          }
        }
      } else {
        if (!quiet) {
          PRINTFB(G, FB_ObjectMap, FB_Details)
          " ObjectMap: Not normalizing: mean = %8.6f and stdev = %8.6f.\n",
              mean, stdev ENDFB(G);
        }
      }
      oMap->Active = true;
    }
  }
  FreeP(point);
  FreeP(sfidx);
  FreeP(atom_sf);
  FreeP(b_factor);
  FreeP(occup);
  return (c);
}
int SelectorMapCoulomb(PyMOLGlobals* G, int sele1, ObjectMapState* oMap,
                       float cutoff, int state, int neutral, int shift,
                       float shift_power) {
  CSelector* I = G->Selector;
  float* v2;
  int a, b, c, j;
  int at;
  int s, idx;
  AtomInfoType* ai;
  ObjectMolecule* obj;
  CoordSet* cs;
  int state1, state2;
  int once_flag;
  int n_at = 0;
  double tot_charge = 0.0;
  float* point = nullptr;
  float* charge = nullptr;
  int n_point = 0;
  int n_occur;
  float *v0, *v1;
  float c_factor = 1.0F;
  float cutoff_to_power = 1.0F;
  const float _1 = 1.0F;

  if (shift) cutoff_to_power = (float)pow(cutoff, shift_power);

  c_factor = SettingGetGlobal_f(G, cSetting_coulomb_units_factor) /
             SettingGetGlobal_f(G, cSetting_coulomb_dielectric);

  c = 0;
  SelectorUpdateTable(G, state, -1);

  point = VLAlloc(float, I->Table.size() * 3);
  charge = VLAlloc(float, I->Table.size());

  /* first count # of times each atom appears */

  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    at = I->Table[a].atom;
    obj = I->Obj[I->Table[a].model];
    s = obj->AtomInfo[at].selEntry;
    ai = obj->AtomInfo + at;
    if (SelectorIsMember(G, s, sele1)) {
      n_occur = 0;
      /* count */
      once_flag = true;
      for (state2 = 0; state2 < obj->NCSet; state2++) {
        if (state < 0) once_flag = false;
        if (!once_flag)
          state1 = state2;
        else
          state1 = state;
        if (state1 < obj->NCSet)
          cs = obj->CSet[state1];
        else
          cs = nullptr;
        if (cs) {
          idx = cs->atmToIdx(at);
          if (idx >= 0) {
            n_occur++;
            n_at++;
          }
        }
        if (once_flag) break;
      }
      /* copy */
      if (n_occur) {
        once_flag = true;
        for (state2 = 0; state2 < obj->NCSet; state2++) {
          if (state < 0) once_flag = false;
          if (!once_flag)
            state1 = state2;
          else
            state1 = state;
          if (state1 < obj->NCSet)
            cs = obj->CSet[state1];
          else
            cs = nullptr;
          if (cs) {
            idx = cs->atmToIdx(at);
            if (idx >= 0) {
              VLACheck(point, float, 3 * n_point + 2);
              VLACheck(charge, float, n_point);
              v0 = cs->coordPtr(idx);
              v1 = point + 3 * n_point;
              copy3f(v0, v1);
              charge[n_point] = ai->partialCharge * ai->q / n_occur;

              tot_charge += charge[n_point];
              n_point++;
            }
          }
          if (once_flag) break;
        }
      }
    }
  }

  PRINTFB(G, FB_Selector, FB_Details)
  " %s: Total charge is %0.3f for %d points (%d atoms).\n", __func__,
      tot_charge, n_point, n_at ENDFB(G);

  if (neutral && (fabs(tot_charge) > R_SMALL4)) {
    float adjust;

    adjust = (float)(-tot_charge / n_point);

    for (a = 0; a < n_point; a++) {
      charge[a] += adjust;
    }

    PRINTFB(G, FB_Selector, FB_Details)
    " %s: Setting net charge to zero...\n", __func__ ENDFB(G);
  }

  for (a = 0; a < n_point; a++) { /* premultiply c_factor by charges */
    charge[a] *= c_factor;
  }

  /* now create and apply voxel map */
  c = 0;
  if (n_point) {
    int* min = oMap->Min;
    int* max = oMap->Max;
    CField* data = oMap->Field->data.get();
    CField* points = oMap->Field->points.get();
    float dist;

    if (cutoff > 0.0F) { /* we are using a cutoff */
      if (shift) {
        PRINTFB(G, FB_Selector, FB_Details)
        " %s: Evaluating local Coulomb potential for grid (shift=%0.2f)...\n",
            __func__, cutoff ENDFB(G);
      } else {
        PRINTFB(G, FB_Selector, FB_Details)
        " %s: Evaluating Coulomb potential for grid (cutoff=%0.2f)...\n",
            __func__, cutoff ENDFB(G);
      }

      std::unique_ptr<MapType> map(
          new MapType(G, -(cutoff), point, n_point, nullptr));
      if (map) {
        float dx, dy, dz;
        float cut = cutoff;
        float cut2 = cutoff * cutoff;

        for (a = min[0]; a <= max[0]; a++) {
          OrthoBusyFast(G, a - min[0], max[0] - min[0] + 1);
          for (b = min[1]; b <= max[1]; b++) {
            for (c = min[2]; c <= max[2]; c++) {
              F3(data, a, b, c) = 0.0F;
              v2 = F4Ptr(points, a, b, c, 0);
              {
                {
                  for (const auto j : MapEIter(*map, v2)) {
                    v1 = point + 3 * j;
                    while (1) {
                      dx = v1[0] - v2[0];
                      dy = v1[1] - v2[1];
                      dx = (float)fabs(dx);
                      dy = (float)fabs(dy);
                      if (dx > cut) break;
                      dz = v1[2] - v2[2];
                      dx = dx * dx;
                      if (dy > cut) break;
                      dz = (float)fabs(dz);
                      dy = dy * dy;
                      if (dz > cut) break;
                      dx = dx + dy;
                      dz = dz * dz;
                      if (dx > cut2) break;
                      dy = dx + dz;
                      if (dy > cut2) break;
                      dist = (float)sqrt1f(dy);

                      if (dist > R_SMALL4) {
                        if (shift) {
                          if (dist < cutoff) {
                            F3(data, a, b, c) +=
                                (charge[j] / dist) *
                                (_1 - (float)pow(dist, shift_power) /
                                          cutoff_to_power);
                          }
                        } else {
                          F3(data, a, b, c) += charge[j] / dist;
                        }
                      }

                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else {
      float* v1;
      PRINTFB(G, FB_Selector, FB_Details)
      " %s: Evaluating Coulomb potential for grid (no cutoff)...\n",
          __func__ ENDFB(G);

      for (a = min[0]; a <= max[0]; a++) {
        OrthoBusyFast(G, a - min[0], max[0] - min[0] + 1);
        for (b = min[1]; b <= max[1]; b++) {
          for (c = min[2]; c <= max[2]; c++) {
            F3(data, a, b, c) = 0.0F;
            v1 = point;
            v2 = F4Ptr(points, a, b, c, 0);
            for (j = 0; j < n_point; j++) {
              dist = (float)diff3f(v1, v2);
              v1 += 3;
              if (dist > R_SMALL4) {
                F3(data, a, b, c) += charge[j] / dist;
              }
            }
          }
        }
      }
    }
    oMap->Active = true;
  }
  VLAFreeP(point);
  VLAFreeP(charge);
  return (1);
}
