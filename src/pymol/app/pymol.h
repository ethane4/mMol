
/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC.
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

#ifndef _H_CPyMOL
#define _H_CPyMOL

/* WARNING: DO NOT DEVELOP APPLICATIONS AGAINST THIS PyMOL_* API.
 *
 * This is a internal interface to PyMOL for the exclusive use of
 * Schrodinger, LLC in development of PyMOL itself, for
 * development of "wrapped" PyMOL products, and for use as foundation
 * layer for support of stable public APIs.
 *
 * Schrodinger, LLC changes this interface regularly and without
 * notice.  It may even vanish altogether.  Code you develop against
 * this interface is likely to be fragile, time-consuming, and
 * expensive to maintain.
 *
 * Our recommended public APIs for PyMOL development are Python-based:
 * (#1) the old "from pymol import cmd" module-based interface, or
 * (#2) the new "from pymol2 import PyMOL; cmd=PyMOL().cmd"
 * instance-based interface.  If you need a C, C++, Java, or
 * ActiveX/COM interface or widget for molecular visualization, then
 * please contact sales@pymol.org about obtaining access to our
 * independent developer and/or OEM product APIs (e.g. JyMOL).
 *
 * If you feel that you absolutely must rely upon this PyMOL_* API,
 * then please be sure to create your own lightweight wrapper layer
 * around it so that the rest of your code will be relatively
 * insensitive to changes made here.  Failure to do so will almost
 * certainly result in code that is impossible to maintain over time.
 *
 * You have been warned!
 */

#include "ov/ov_returns.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PYMOL_BUTTON_DOWN 0
#define PYMOL_BUTTON_UP 1

#define PYMOL_KEY_F1 1
#define PYMOL_KEY_F2 2
#define PYMOL_KEY_F3 3
#define PYMOL_KEY_F4 4
#define PYMOL_KEY_F5 5
#define PYMOL_KEY_F6 6
#define PYMOL_KEY_F7 7
#define PYMOL_KEY_F8 8
#define PYMOL_KEY_F9 9
#define PYMOL_KEY_F10 10
#define PYMOL_KEY_F11 11
#define PYMOL_KEY_F12 12
#define PYMOL_KEY_LEFT 100
#define PYMOL_KEY_UP 101
#define PYMOL_KEY_RIGHT 102
#define PYMOL_KEY_DOWN 103
#define PYMOL_KEY_PAGE_UP 104
#define PYMOL_KEY_PAGE_DOWN 105
#define PYMOL_KEY_HOME 106
#define PYMOL_KEY_END 107
#define PYMOL_KEY_INSERT 108

#define PYMOL_BUTTON_LEFT 0
#define PYMOL_BUTTON_MIDDLE 1
#define PYMOL_BUTTON_RIGHT 2
#define PYMOL_BUTTON_SCROLL_FORWARD 3
#define PYMOL_BUTTON_SCROLL_REVERSE 4

#define PYMOL_MODIFIER_SHIFT 1
#define PYMOL_MODIFIER_CTRL 2
#define PYMOL_MODIFIER_ALT 4

#define PYMOL_FALSE 0
#define PYMOL_TRUE 1

#define PYMOL_DEFAULT -1

#define PYMOL_PROGRESS_SLOW 0
#define PYMOL_PROGRESS_MED 2
#define PYMOL_PROGRESS_FAST 4

#define PYMOL_PROGRESS_SIZE 6
#define PYMOL_RESHAPE_SIZE 5

/* configuration */

struct CPyMOLOptions;

CPyMOLOptions* PyMOLOptions_New(void);
void PyMOLOptions_Free(CPyMOLOptions* option);

/* PyMOL instance type */

struct CPyMOL;

#define ATOM_PROP_MODEL 0
#define ATOM_PROP_INDEX 1
#define ATOM_PROP_TYPE 2
#define ATOM_PROP_NAME 3
#define ATOM_PROP_RESN 4
#define ATOM_PROP_RESI 5
#define ATOM_PROP_RESV 6
#define ATOM_PROP_CHAIN 7
#define ATOM_PROP_ALT 8
#define ATOM_PROP_SEGI 9
#define ATOM_PROP_ELEM 10
#define ATOM_PROP_SS 11
#define ATOM_PROP_TEXT_TYPE 12
#define ATOM_PROP_CUSTOM 13
#define ATOM_PROP_LABEL 14
#define ATOM_PROP_NUMERIC_TYPE 15
#define ATOM_PROP_Q 16
#define ATOM_PROP_B 17
#define ATOM_PROP_VDW 18
#define ATOM_PROP_ELEC_RADIUS 19
#define ATOM_PROP_PARTIAL_CHARGE 20
#define ATOM_PROP_FORMAL_CHARGE 21
#define ATOM_PROP_STEREO 22
#define ATOM_PROP_CARTOON 23
#define ATOM_PROP_COLOR 24
#define ATOM_PROP_ID 25
#define ATOM_PROP_RANK 26
#define ATOM_PROP_FLAGS 27
#define ATOM_PROP_GEOM 28
#define ATOM_PROP_VALENCE 29
#define ATOM_PROP_X 30
#define ATOM_PROP_Y 31
#define ATOM_PROP_Z 32
#define ATOM_PROP_SETTINGS 33
#define ATOM_PROP_PROPERTIES 34
#define ATOM_PROP_ONELETTER 40
#define ATOM_PROP_EXPLICIT_DEGREE 41
#define ATOM_PROP_EXPLICIT_VALENCE 42

/* return status values */

#define PyMOLstatus_YES 1
#define PyMOLstatus_NO 0

#define PyMOLstatus_SUCCESS 0
#define PyMOLstatus_FAILURE -1

/* return types */

typedef int PyMOLstatus;

typedef struct {
  PyMOLstatus status;
} PyMOLreturn_status;

typedef struct {
  PyMOLstatus status;
  int value;
} PyMOLreturn_int;

typedef struct {
  PyMOLstatus status;
  float value;
} PyMOLreturn_float;

typedef struct {
  PyMOLstatus status;
  int size;
  float* array;
} PyMOLreturn_float_array;

typedef struct {
  PyMOLstatus status;
  int size;
  int* array;
} PyMOLreturn_int_array;

typedef struct {
  PyMOLstatus status;
  int size;
  char** array;
} PyMOLreturn_string_array;

#define PYMOL_RETURN_VALUE_IS_STRING 0x01
#define PYMOL_RETURN_VALUE_IS_INT 0x02
#define PYMOL_RETURN_VALUE_IS_FLOAT 0x04
#define PYMOL_RETURN_VALUE_IS_FLOAT_ARRAY 0x08

typedef struct {
  PyMOLstatus status;
  short int type;
  char* string;
  int int_value;
  float float_value;
  int array_length;
  float* float_array;
} PyMOLreturn_value;

typedef void PyMOLModalDrawFn(void* G);

/* creation */

CPyMOL* PyMOL_NewWithOptions(const CPyMOLOptions* option);

/* destruction */

void PyMOL_Free(CPyMOL* I);

/* starting and stopping */

void PyMOL_Start(CPyMOL* I);
void PyMOL_StartWithPython(CPyMOL* I);
void PyMOL_SetPythonInitStage(CPyMOL* I, int stage);
void PyMOL_Stop(CPyMOL* I);

/* upstream invalidation and configuration events */

void PyMOL_NeedFakeDrag(CPyMOL* I);
void PyMOL_NeedRedisplay(CPyMOL* I);
void PyMOL_NeedSwap(CPyMOL* I);
void PyMOL_SetClickReady(CPyMOL* I, const char* name, int index, int button,
                         int mod, int x, int y, const float* pos, int state,
                         int bond = -1);
void PyMOL_SetPassive(CPyMOL* I, int onOff);
void PyMOL_NeedReshape(CPyMOL* I, int mode, int x, int y, int width,
                       int height);

/* valid context management */

void PyMOL_PushValidContext(CPyMOL* I);
void PyMOL_PopValidContext(CPyMOL* I);

/* methods requiring a valid OpenGL context*/

void PyMOL_Draw(CPyMOL* I);

/* methods that do not require a valid OpenGL context */

void PyMOL_Key(CPyMOL* I, unsigned char k, int x, int y, int modifiers);
void PyMOL_Special(CPyMOL* I, int k, int x, int y, int modifiers);
void PyMOL_Reshape(CPyMOL* I, int width, int height, int force);
void PyMOL_Drag(CPyMOL* I, int x, int y, int modifiers);
void PyMOL_Button(CPyMOL* I, int button, int state, int x, int y,
                  int modifiers);
int PyMOL_Idle(CPyMOL* I); /* return true if PyMOL is busy doing real
                              work (not simply racing the CPU) */

void PyMOL_ExpireIfIdle(CPyMOL* I); /* auto-termination for command-line mode */

typedef void PyMOLSwapBuffersFn(void);

void PyMOL_SwapBuffers(
    CPyMOL* I); /* only works if above function has been set */
void PyMOL_SetDefaultMouse(CPyMOL* I);

/* host query methods */

int PyMOL_GetRedisplay(CPyMOL* I, int reset);

int PyMOL_GetIdleAndReady(CPyMOL* I);

/* int array results */

/* string results */

/* asynchronous processing (only useful for Python-based multithreaded builds
 * right now) */

int PyMOL_GetBusy(CPyMOL* I, int reset);
void PyMOL_SetBusy(CPyMOL* I, int value);

void PyMOL_ResetProgress(CPyMOL* I);

void PyMOL_SetProgress(CPyMOL* I, int offset, int current, int range);

int PyMOL_GetProgress(CPyMOL* I, int* progress, int reset);

int PyMOL_GetInterrupt(CPyMOL* I, int reset);
void PyMOL_SetInterrupt(CPyMOL* I, int value);

/* modal updates -- PyMOL is busy with some complex task, but we have
   to return control to the host in order to get a valid draw callback */

int PyMOL_GetModalDraw(CPyMOL* I);
void PyMOL_SetModalDraw(CPyMOL* I,
                        PyMOLModalDrawFn* fn); /* for internal use only */

/* developer/transient privates */

struct PyMOLGlobals* PyMOL_GetGlobals(CPyMOL* I);
struct PyMOLGlobals** PyMOL_GetGlobalsHandle(CPyMOL* I);

/* for python integration */

/* for Jmol integration */

/* Command API */

PyMOLreturn_status PyMOL_CmdLoad(CPyMOL* I, const char* content,
                                 const char* content_type,
                                 const char* content_format,
                                 const char* object_name, int state,
                                 int discrete, int finish, int quiet,
                                 int multiplex, int zoom);

/* releasing returned values */

PyMOLreturn_status PyMOL_ZoomScene(CPyMOL* I, float scale);
PyMOLreturn_status PyMOL_TranslateScene(CPyMOL* I, float x, float y, float z);

typedef struct _P_AtomProperty {
  int id;
  short Ptype;
  int offset;
  int maxlen;
} AtomPropertyInfo;

AtomPropertyInfo* PyMOL_GetAtomPropertyInfo(CPyMOL* I,
                                            const char* atompropname);

#ifdef __cplusplus
}
#endif
#endif
