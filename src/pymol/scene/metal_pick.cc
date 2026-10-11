/*
 * Screen-space atom picking for the Metal backend -- see MetalPick.h.
 */

#include "pymol/scene/metal_pick.h"

#include <algorithm>
#include <cfloat>

#include "pymol/objects/atom_info.h"
#include "pymol/objects/coord_set.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/reps/rep.h"
#include "pymol/settings/setting.h"
#include "pymol/geometry/vector.h"

namespace {

/// An atom that projected inside the pick radius.
struct PickCand {
  float d2;
  float
      w2;  // cursor-to-atom distance squared, in Angstroms at the atom's depth
  float depth;
  ObjectMolecule* obj;
  int atm;
  float sx, sy;
};

/**
 * True if `ai` has geometry a pick can land on.
 *
 * Cartoon and ribbon are special: showing them ORs their bit onto EVERY atom of
 * the object (side chains and solvent included), but the spline only runs
 * through the guide atoms, so only those are hittable. Mirrors the _DRAWN_REPS
 * selection in modules/mmol/pick.py.
 */
inline bool is_drawn(const AtomInfoType* ai, int rep_mask, int guide_rep_mask) {
  if (ai->visRep & rep_mask) return true;
  return (ai->visRep & guide_rep_mask) && (ai->flags & cAtomFlag_guide);
}

}  // namespace

MetalPickHit MetalPickAtom(PyMOLGlobals* G,
                           const std::vector<ObjectMolecule*>& objects,
                           int state, const MetalPickCamera& cam, float ndc_x,
                           float ndc_y, float max_ndc2, float cluster_ndc2,
                           int rep_mask, int guide_rep_mask) {
  MetalPickHit hit;

  const float r00 = cam.rot[0], r01 = cam.rot[1], r02 = cam.rot[2];
  const float r10 = cam.rot[3], r11 = cam.rot[4], r12 = cam.rot[5];
  const float r20 = cam.rot[6], r21 = cam.rot[7], r22 = cam.rot[8];
  const float tx = cam.pos[0], ty = cam.pos[1], tz = cam.pos[2];
  const float ox = cam.origin[0], oy = cam.origin[1], oz = cam.origin[2];
  // A degenerate slab (front >= back) means the view layout didn't carry usable
  // clip planes; culling against it would make nothing pickable at all.
  const bool clipped = cam.clip_back > cam.clip_front;

  if (!(cam.tan_half > 0.f) || !(cam.aspect > 0.f)) return hit;

  std::vector<PickCand> cands;

  constexpr float kReachA = 2.0f;
  for (auto* obj : objects) {
    if (!obj) continue;
    // all_states draws every state at once: pick from all of them, or a hover
    // over state 3 snaps to the nearest atom of the current state instead.
    const bool every = state == -2 && ObjectGetCurrentState(obj, false) == -1;
    for (int st = every ? 0 : state, st_end = every ? obj->NCSet : state + 1;
         st < st_end; ++st) {
      const CoordSet* cs = obj->getCoordSet(st);
      if (!cs) continue;

      // Same transform CoordSetGetAtomTxfVertex applies, hoisted out of the
      // loop: the pick has to agree with where the renderer put the atom, and
      // for a moved object (non-identity TTT) that is nowhere near the raw
      // coordinate.
      const double* mat = nullptr;
      if (!cs->Matrix.empty() && SettingGet<int>(*cs, cSetting_matrix_mode) > 0)
        mat = cs->Matrix.data();
      const bool use_ttt = obj->TTTFlag != 0;

      const AtomInfoType* atom_info = obj->AtomInfo.data();
      const int n_index = cs->getNIndex();

      for (int idx = 0; idx < n_index; ++idx) {
        const int atm = cs->IdxToAtm[idx];
        const AtomInfoType* ai = atom_info + atm;
        if (!is_drawn(ai, rep_mask, guide_rep_mask)) continue;

        float v[3];
        copy3f(cs->coordPtr(idx), v);
        if (mat) transform44d3f(mat, v, v);
        if (use_ttt) transformTTT44f3f(obj->TTT, v, v);

        const float dx = v[0] - ox;
        const float dy = v[1] - oy;
        const float dz = v[2] - oz;

        // eye = R*(model-origin) + pos; the camera looks down -Z.
        const float depth = -(r20 * dx + r21 * dy + r22 * dz + tz);
        if (depth <= 0.01f) continue;
        // An atom clipped away isn't visible, so it must not be selectable.
        if (clipped && (depth < cam.clip_front || depth > cam.clip_back))
          continue;

        const float half_h =
            cam.ortho_half > 0.f ? cam.ortho_half : depth * cam.tan_half;
        const float sx = (r00 * dx + r01 * dy + r02 * dz + tx) /
                         (half_h * cam.aspect);  // NDC x, +1 = right
        const float sy =
            (r10 * dx + r11 * dy + r12 * dz + ty) / half_h;  // NDC y, +1 = up

        const float ex = sx - ndc_x;
        const float ey = sy - ndc_y;
        const float d2 = ex * ex + ey * ey;
        const float wx = ex * half_h * cam.aspect, wy = ey * half_h;
        const float w2 = wx * wx + wy * wy;
        // within the screen radius, or within kReachA of the atom when zoomed
        // in
        if (d2 > max_ndc2 && w2 > kReachA * kReachA) continue;

        ++hit.ncand;
        cands.push_back({d2, w2, depth, obj, atm, sx, sy});
      }
    }
  }

  if (cands.empty()) return hit;

  // Of the atoms whose footprint (kCoverA around the atom, at its depth) the
  // cursor is on, the front-most is the one visible there. A fixed screen-space
  // cluster instead let a nearby front atom win -- in a helix, residue i+4 over
  // the residue under the cursor -- and more so the further out the zoom.
  // With none under the cursor, the nearest on screen.
  constexpr float kCoverA = 1.0f;
  std::stable_sort(
      cands.begin(), cands.end(),
      [](const PickCand& a, const PickCand& b) { return a.d2 < b.d2; });
  const PickCand* best = nullptr;
  for (const auto& c : cands)
    if (c.w2 < kCoverA * kCoverA && (!best || c.depth < best->depth)) best = &c;
  if (!best) best = &cands.front();
  (void)cluster_ndc2;

  hit.obj = best->obj;
  hit.atm = best->atm;
  hit.d2 = best->d2;
  hit.sx = best->sx;
  hit.sy = best->sy;
  return hit;
}
