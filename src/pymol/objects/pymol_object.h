
/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
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
#ifndef _H_PyMOLObject
#define _H_PyMOLObject

/* literally a 3-D object...also an object object */

#include "pymol/base/base.h"
#include "pymol/reps/rep.h"
#include "pymol/settings/setting.h"
#include "pymol/symmetry/symmetry.h"
#include "pymol/base/pymol_globals.h"
#include "pymol/scene/view.h"
#include "pymol/base/word.h"
#include "pymol/math_defines.h"
#include "pymol/base/vla.h"
#include "pymol/graphics/render_context.h"
#include "pymol/scene/mview_action.h"

#include <string>
#include <type_traits>

using ObjectNameType = char[WordLength];

enum cObject_t : int {
  cObjectMolecule = 1,
  cObjectMap = 2,
  cObjectMesh = 3,
  cObjectMeasurement = 4,
  cObjectCallback = 5,
  cObjectCGO = 6,
  cObjectSurface = 7,
  cObjectGadget = 8,
  cObjectCalculator = 9,
  cObjectSlice = 10,
  cObjectAlignment = 11,
  cObjectGroup = 12,
  cObjectVolume = 13,
  cObjectCurve = 14,
  // Gizmo is UI object in world coordinate space (in 3D)
  cObjectGizmo = 15,
};

/*
   the object base class is in the process of being converted to support
   states explicitly (an unfortunate early omission), which will allow
   for simplified implementation of future multi-state objects.
 */

struct CObjectState {
  PyMOLGlobals* G;
  std::vector<double> Matrix;
  std::vector<double> InvMatrix;
  CObjectState() : G(nullptr) {};
  CObjectState(PyMOLGlobals* G) : G(G) {};
};

namespace pymol {
struct CObject {
  PyMOLGlobals* G = nullptr;
  cObject_t type;
  ObjectNameType Name{};
  int Color = 0;
  int visRep = 0;
  float ExtentMin[3]{}, ExtentMax[3]{};
  int ExtentFlag = false, TTTFlag = false;
  float TTT[16]{}; /* translate, transform, translate matrix (to apply when
                      rendering) */
  pymol::copyable_ptr<CSetting> Setting;
  int Enabled = 0;                /* read-only... maintained by Scene */
  pymol::vla<CViewElem> ViewElem; /* for animating objects via the TTT */

  /* not pickled */
  int grid_slot = 0;
  CGO* gridSlotSelIndicatorsCGO = nullptr;
  int Grabbed = 0;
  // Movie frame at which the ViewElem per-frame state track was last written
  // into cSetting_state (see ObjectPrepareContext). Re-applying it on EVERY
  // render at a STATIC frame would clobber an explicit per-object `set state`
  // (the Objects-panel model inspection, which is decoupled from the movie) —
  // freezing the viewport on one model while a movie sits idle (issue #132).
  // The movie's own play / scrub / export all ADVANCE the frame, so the sweep
  // still drives them. -1 = not yet applied.
  int LastStateFrame = -1;

  // methods
  StateIndex_t getCurrentState() const;
  CObjectState* getObjectState(StateIndex_t state);
  const CObjectState* getObjectState(StateIndex_t state) const {
    return const_cast<CObject*>(this)->getObjectState(state);
  }

 protected:
  CObject(PyMOLGlobals* G);

  /// @pre `0 <= state && state < getNFrame()`
  /// @return nullptr if state is empty (not active)
  virtual CObjectState* _getObjectState(int state) { return nullptr; }

 public:
  virtual ~CObject();

  /**
   * Gets a copy TTT matrix of PyMOL Object, disregarding TTT flag
   * @return copy of TTT matrix
   * Note: See also ObjectGetTTT (by-pointer retrieval)
   */
  pymol::Mat4 getTTT() const {
    pymol::Mat4 ret_ttt;
    std::copy_n(std::begin(TTT), 16, ret_ttt.begin());
    return ret_ttt;
  }

  /**
   * Method version of ObjectSetName
   * @param name new object Name
   */
  void setName(pymol::zstring_view name);

  virtual void update() {}
  virtual void render(RenderInfo* info);
  virtual void invalidate(cRep_t rep, cRepInv_t level, int state) {}
  virtual int getNFrame() const { return 1; }
  virtual std::string describeElement(int index) const;
  virtual char* getCaption(char* ch, int len) const { return nullptr; };
  virtual pymol::copyable_ptr<CSetting>* getSettingHandle(int state);
  virtual CObject* clone() const { return nullptr; };
  virtual CSymmetry const* getSymmetry(int state = 0) const { return nullptr; }
  virtual bool setSymmetry(CSymmetry const&, int state = 0) { return false; }
  virtual pymol::RenderContext getRenderContext() const {
    return RenderContext::Camera;
  }
};
}  // namespace pymol
int ObjectCopyHeader(pymol::CObject* I, const pymol::CObject* src);
void ObjectSetName(pymol::CObject* I, const char* name);
bool ObjectMakeValidName(char* name);
void ObjectMakeValidName(PyMOLGlobals* G, char* name, bool quiet = false);
void ObjectPurgeSettings(pymol::CObject* I);
void ObjectUseColor(pymol::CObject* I);
void ObjectUseColorCGO(CGO* cgo, pymol::CObject* I);
void ObjectSetRepVisMask(pymol::CObject* I, int repmask, int value);
void ObjectPrepareContext(pymol::CObject* I, RenderInfo* info);
void ObjectSetTTT(pymol::CObject* I, const float* ttt, int state, int store);
int ObjectGetTTT(pymol::CObject* I, const float** ttt, int state);
int ObjectGetTotalMatrix(pymol::CObject* I, int state, int history,
                         double* matrix);
void ObjectCombineTTT(pymol::CObject* I, const float* ttt, int reverse_order,
                      int store);
void ObjectTranslateTTT(pymol::CObject* T, const float* v, int store);
void ObjectSetTTTOrigin(pymol::CObject* I, float* origin);
void ObjectResetTTT(pymol::CObject* I, int store);
PyObject* ObjectAsPyList(const pymol::CObject* I);
int ObjectFromPyList(PyMOLGlobals* G, PyObject* list, pymol::CObject* I);
int ObjectGetCurrentState(const pymol::CObject* I, int ignore_all_states);
void ObjectAdjustStateRebuildRange(pymol::CObject* I, int* start, int* stop);
int ObjectMotion(pymol::CObject* I, MViewAction action, int first, int last,
                 float power, float bias, int simple, float linear, int wrap,
                 int hand, int window, int cycles, int state, int quiet);
int ObjectGetSpecLevel(pymol::CObject* I, int frame);
void ObjectMotionTrim(pymol::CObject* I, int n_frame);
void ObjectDrawViewElem(pymol::CObject* I, BlockRect* rect, int frames,
                        CGO* orthoCGO);
void ObjectStateInit(PyMOLGlobals* G, CObjectState* I);
void ObjectStatePurge(CObjectState* I);
int ObjectStateSetMatrix(CObjectState* I, const double* matrix);
double* ObjectStateGetMatrix(CObjectState* I);
const double* ObjectStateGetMatrix(const CObjectState* I);
const double* ObjectStateGetInvMatrix(const CObjectState* I);
void ObjectStateTransformMatrix(CObjectState* I, const double* matrix);
void ObjectStateResetMatrix(CObjectState* I);
PyObject* ObjectStateAsPyList(CObjectState* I);
int ObjectStateFromPyList(PyMOLGlobals* G, PyObject* list, CObjectState* I);
int ObjectStatePushAndApplyMatrix(CObjectState* I, RenderInfo* info);
void ObjectStatePopMatrix(CObjectState* I, RenderInfo* info);
void ObjectStateLeftCombineMatrixR44d(CObjectState* I, const double* matrix);
int ObjectMotionModify(pymol::CObject* I, ViewElemAction action, int index,
                       int count, int target, int freeze, int localize);
void ObjectMotionReinterpolate(pymol::CObject* I);
int ObjectMotionGetLength(pymol::CObject* I);

struct CObjectUpdateThreadInfo;

#define cObjectTypeAll 0
#define cObjectTypeObjects 1
#define cObjectTypeSelections 2
#define cObjectTypePublic 3
#define cObjectTypePublicObjects 4
#define cObjectTypePublicSelections 5
#define cObjectTypePublicNonGroupObjects 6
#define cObjectTypePublicGroupObjects 7
#define cObjectTypeNonGroupObjects 8
#define cObjectTypeGroupObjects 9
/* Note: public objects are ones that do not start with "_" */

// object and object-state level setting
template <typename V>
void SettingSet(int index, V value, pymol::CObject* obj,
                StateIndex_t state = cStateAll) {
  auto handle = obj->getSettingHandle(state);
  if (handle) SettingSet(obj->G, *handle, index, value);
}

template <typename V>
V SettingGet(const pymol::CObject& obj, int index) {
  using T = typename std::conditional<std::is_enum<V>::value, int, V>::type;
  return static_cast<V>(
      SettingGet<T>(obj.G, obj.Setting.get(), nullptr, index));
}

/**
 * Similar to StateIterator but easily useable
 * with ranged-for loops
 */
class StateIteratorV2 {
  StateIndex_t m_beg;
  StateIndex_t m_end;
  StateIndex_t m_cur;

 public:
  StateIteratorV2(pymol::CObject*, StateIndex_t);

  /**
   * @return begin state iterator
   */
  StateIndex_t* begin();

  /**
   * @return end state iterator
   */
  StateIndex_t* end();
};

#endif
