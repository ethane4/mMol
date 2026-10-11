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

// Reading AMBER topology, trajectory and restart files.

#include "pymol/base/parse.h"
#include "pymol/base/err.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/executive/executive.h"
#include "pymol/base/file.h"
#include "pymol/base/feedback.h"
#include "pymol/objects/object_molecule_internal.h"

static const char* skip_fortran(int num, int per_line, const char* p) {
  int a, b;
  b = 0;
  for (a = 0; a < num; a++) {
    if ((++b) == per_line) {
      b = 0;
      p = nextline(p);
    }
  }
  if (b || (!num)) p = nextline(p);
  return (p);
}

/**
 * Selection mapping for `load_traj` to only load coordinates for a subset of
 * atoms. Returns a mapping of old to new coordinate indices and modifies the
 * index-related members of `cs`.
 *
 * @param[in] obj Object Molecule for evaluation of the selection
 * @param[in,out] cs Candidate coordinate set (with valid index arrays)
 * @param[in] selection Named selection
 * @return index mapping, or nullptr if `selection` is not valid
 */
std::unique_ptr<int[]> LoadTrajSeleHelper(const ObjectMolecule* obj,
                                          CoordSet* cs, const char* selection) {
  auto G = obj->G;
  int sele0 = SelectorIndexByName(G, selection);

  // We can skip "all" (0) here since it's a trivial case
  if (sele0 <= 0) {
    return nullptr;
  }

  auto xref = std::unique_ptr<int[]>(new int[cs->NIndex]);

  int idx_new = 0;

  for (int idx = 0; idx < cs->NIndex; ++idx) {
    auto atm = cs->IdxToAtm[idx];
    if (SelectorIsMember(G, obj->AtomInfo[atm].selEntry, sele0)) {
      cs->IdxToAtm[idx_new] = atm;
      cs->AtmToIdx[atm] = idx_new;
      xref[idx] = idx_new;
      ++idx_new;
    } else {
      cs->AtmToIdx[atm] = -1;
      xref[idx] = -1;
    }
  }

  cs->NIndex = idx_new;
  cs->IdxToAtm.resize(cs->NIndex);
  cs->Coord.resize(cs->NIndex * 3);

  return xref;
}
ObjectMolecule* ObjectMoleculeLoadTRJFile(PyMOLGlobals* G, ObjectMolecule* I,
                                          const char* fname, int frame,
                                          int interval, int average, int start,
                                          int stop, int max, const char* sele,
                                          int image, const float* shift,
                                          int quiet) {
  FILE* f;
  char* buffer;
  const char* p;
  char cc[MAXLINELEN];
  int n_read;
  int to_go;
  int skip_first_line = true;
  int periodic = false;
  int angles = true;
  float f0, f1, f2, *fp;
  float box[3], angle[3];
  float r_cent[3], r_trans[3];
  int r_act, r_val, r_cnt;
  float *r_fp_start = nullptr, *r_fp_stop = nullptr;
  int a, b, c, i;
  int zoom_flag = false;
  int cnt = 0;
  int n_avg = 0;
  int icnt;
  int ncnt = 0;
  float zerovector[3] = {0.0, 0.0, 0.0};
  CoordSet* cs = nullptr;

  if (I->DiscreteFlag) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " %s: Discrete objects not supported\n", __func__ ENDFB(G);
    return I;
  }

  if (!shift) shift = zerovector;
  if (interval < 1) interval = 1;

  icnt = interval;
#define BUFSIZE 4194304
#define GETTING_LOW 10000

  f = pymol_fopen(fname, "rb");
  if (!f) {
    ErrMessage(G, __func__, "Unable to open file!");
  } else {
    if (I->CSTmpl) {
      cs = CoordSetCopy(I->CSTmpl);
    } else if (I->NCSet > 0) {
      cs = CoordSetCopy(I->CSet[0]);
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjMolLoadTRJFile: Missing topology" ENDFB(G);
      return (I);
    }

    auto xref = LoadTrajSeleHelper(I, cs, sele);

    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    " ObjMolLoadTRJFile: Loading from \"%s\".\n", fname ENDFB(G);
    buffer = pymol::malloc<char>(BUFSIZE + 1); /* 1 MB read buffer */
    p = buffer;
    buffer[0] = 0;
    n_read = 0;
    to_go = 0;
    a = 0;
    b = 0;
    c = 0;
    f1 = 0.0;
    f2 = 0.0;
    while (1) {
      to_go = n_read - (p - buffer);
      if (to_go < GETTING_LOW)
        if (!feof(f)) {
          if (to_go) memcpy(buffer, p, to_go);
          n_read = fread(buffer + to_go, 1, BUFSIZE - to_go, f);
          n_read = to_go + n_read;
          buffer[n_read] = 0;
          p = buffer;
          if (skip_first_line) {
            p = nextline(p);
            skip_first_line = false;
          }
          to_go = n_read - (p - buffer);
        }
      if (!to_go) break;
      p = ncopy(cc, p, 8);
      if ((++b) == 10) {
        b = 0;
        p = nextline(p);
      }
      f0 = f1;
      f1 = f2;
      if (sscanf(cc, "%f", &f2) == 1) {
        if ((++c) == 3) {
          c = 0;
          if ((cnt + 1) >= start) {
            if (icnt <= 1) {
              if (xref) {
                if (xref[a] >= 0)
                  fp = cs->coordPtr(xref[a]);
                else
                  fp = nullptr;
              } else {
                fp = cs->coordPtr(a);
              }
              if (fp) {
                if (n_avg) {
                  *(fp++) += f0;
                  *(fp++) += f1;
                  *(fp++) += f2;
                } else {
                  *(fp++) = f0;
                  *(fp++) = f1;
                  *(fp++) = f2;
                }
              }
            }
          }
          if ((++a) == I->NAtom) {
            cnt++;
            a = 0;
            if (b) p = nextline(p);
            b = 0;

            if (cs->PeriodicBoxType != CoordSet::NoPeriodicity) {
              /* read periodic box */

              c = 0;
              periodic = true;
              angles = true;

              p = ncopy(cc, p, 8);
              if (sscanf(cc, "%f", &box[0]) != 1) periodic = false;
              p = ncopy(cc, p, 8);
              if (sscanf(cc, "%f", &box[1]) != 1) periodic = false;
              p = ncopy(cc, p, 8);
              if (sscanf(cc, "%f", &box[2]) != 1) periodic = false;

              p = ncopy(cc, p, 8);
              if (sscanf(cc, "%f", &angle[0]) != 1) angles = false;

              p = ncopy(cc, p, 8);
              if (sscanf(cc, "%f", &angle[1]) != 1) angles = false;

              p = ncopy(cc, p, 8);
              if (sscanf(cc, "%f", &angle[2]) != 1) angles = false;
              if (periodic) {
                cs->Symmetry = std::make_unique<CSymmetry>(G);
                cs->Symmetry->Crystal.setDims(box);
                if (angles) {
                  cs->Symmetry->Crystal.setAngles(angle);
                }
                p = nextline(p);
                b = 0;
              }

              if (cs->PeriodicBoxType == CoordSet::Octahedral)
                periodic = false; /* can't handle this yet... */
            }

            if ((stop > 0) && (cnt >= stop)) break;
            if (cnt >= start) {
              icnt--;
              if (icnt > 0) {
                PRINTFB(G, FB_ObjectMolecule, FB_Details)
                " ObjectMolecule: skipping set %d...\n", cnt ENDFB(G);
              } else {
                icnt = interval;
                n_avg++;
              }

              if (icnt == interval) {
                if (n_avg < average) {
                  PRINTFB(G, FB_ObjectMolecule, FB_Details)
                  " ObjectMolecule: averaging set %d...\n", cnt ENDFB(G);
                } else {
                  /* compute average */
                  if (n_avg > 1) {
                    fp = cs->Coord.data();
                    for (i = 0; i < cs->NIndex; i++) {
                      *(fp++) /= n_avg;
                      *(fp++) /= n_avg;
                      *(fp++) /= n_avg;
                    }
                  }
                  if (periodic && image) {
                    /* Perform residue-based period image transformation */
                    i = 0;
                    r_cnt = 0;
                    r_act = 0; /* 0 unspec, 1=load, 2=image, 3=leave */
                    r_val = -1;
                    while (r_act != 3) {
                      if (i >= cs->NIndex) {
                        if (r_cnt)
                          r_act = 2;
                        else
                          r_act = 3;
                      }
                      if (r_act == 0) {
                        /* start new residue */
                        r_cnt = 0;
                        r_act = 1; /* now load */
                      }
                      if (r_act == 1) {
                        if (i < cs->NIndex) {
                          /* is there a coordinate for atom? */
                          if (xref) {
                            if (xref[i] >= 0)
                              fp = cs->coordPtr(xref[i]);
                            else
                              fp = nullptr;
                          } else {
                            fp = cs->coordPtr(i);
                          }
                          if (fp) { /* yes there is... */
                            if (r_cnt) {
                              if (r_val != I->AtomInfo[cs->IdxToAtm[i]].resv) {
                                r_act = 2; /* end of residue-> time to image */
                              } else {
                                r_cnt++;
                                r_cent[0] += *(fp++);
                                r_cent[1] += *(fp++);
                                r_cent[2] += *(fp++);
                                r_fp_stop = fp; /* stop here */
                                i++;
                              }
                            } else {
                              r_val = I->AtomInfo[cs->IdxToAtm[i]].resv;
                              r_cnt++;
                              r_fp_start = fp; /* start here */
                              r_cent[0] = *(fp++);
                              r_cent[1] = *(fp++);
                              r_cent[2] = *(fp++);
                              r_fp_stop = fp; /* stop here */
                              i++;
                            }
                          } else {
                            i++;
                          }
                        } else {
                          r_act = 2; /* image */
                        }
                      }

                      if (r_act == 2) { /* time to image */
                        if (r_cnt) {
                          r_cent[0] /= r_cnt;
                          r_cent[1] /= r_cnt;
                          r_cent[2] /= r_cnt;
                          const auto& periodicbox = cs->getSymmetry()->Crystal;
                          transform33f3f(periodicbox.realToFrac(), r_cent,
                                         r_cent);
                          r_trans[0] =
                              fmodf(1000.0F + shift[0] + r_cent[0], 1.0F);
                          r_trans[1] =
                              fmodf(1000.0F + shift[1] + r_cent[1], 1.0F);
                          r_trans[2] =
                              fmodf(1000.0F + shift[2] + r_cent[2], 1.0F);
                          r_trans[0] -= r_cent[0];
                          r_trans[1] -= r_cent[1];
                          r_trans[2] -= r_cent[2];
                          transform33f3f(periodicbox.fracToReal(), r_trans,
                                         r_trans);
                          fp = r_fp_start;
                          while (fp < r_fp_stop) {
                            *(fp++) += r_trans[0];
                            *(fp++) += r_trans[1];
                            *(fp++) += r_trans[2];
                          }
                        }
                        r_act = 0; /* reset */
                        r_cnt = 0;
                      }
                    }
                  }

                  /* add new coord set */
                  cs->invalidateRep(cRepAll, cRepInvRep);
                  if (frame < 0) frame = I->NCSet;
                  if (!I->NCSet) {
                    zoom_flag = true;
                  }

                  VLACheck(I->CSet, CoordSet*, frame);
                  if (I->NCSet <= frame) I->NCSet = frame + 1;
                  delete I->CSet[frame];
                  I->CSet[frame] = cs;
                  ncnt++;

                  if (average < 2) {
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                    " ObjectMolecule: read set %d into state %d...\n", cnt,
                        frame + 1 ENDFB(G);
                  } else {
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                    " ObjectMolecule: averaging set %d...\n", cnt ENDFB(G);
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                    " ObjectMolecule: average loaded into state %d...\n",
                        frame + 1 ENDFB(G);
                  }
                  frame++;
                  cs = CoordSetCopy(cs);
                  n_avg = 0;
                  if ((stop > 0) && (cnt >= stop)) break;
                  if ((max > 0) && (ncnt >= max)) break;
                }
              }
            } else {
              PRINTFB(G, FB_ObjectMolecule, FB_Details)
              " ObjectMolecule: skipping set %d...\n", cnt ENDFB(G);
            }
          }
        }
      } else {
        PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " ObjMolLoadTRJFile-Error: Failed to read expected coordinate value.\n "
        " This traj. does not match the loaded parameter/topology file.\n  "
        "Likely cause: either the atom count or the periodic box settings\n  "
        "are inconsistent between the two files.\n" ENDFB(G);
        break;
      }
    }
    mfree(buffer);
  }
  delete cs;
  SceneChanged(G);
  SceneCountFrames(G);
  if (zoom_flag)
    if (SettingGetGlobal_i(G, cSetting_auto_zoom)) {
      ExecutiveWindowZoom(G, I->Name, 0.0, -1, 0, 0,
                          quiet); /* auto zoom (all states) */
    }

  return (I);
}

ObjectMolecule* ObjectMoleculeLoadRSTFile(PyMOLGlobals* G, ObjectMolecule* I,
                                          const char* fname, int frame,
                                          int quiet, char mode) {
  /*
   * mode = 0: AMBER coordinate/restart file (one frame only)
   * mode = 1: AMBER trajectory
   * mode = 2: AMBER trajectory with box
   *           http://ambermd.org/formats.html
   */
  int ok = true;
  char *buffer, *p;
  char cc[MAXLINELEN];
  float f0, f1, f2, *fp;
  int a, b, c;
  int zoom_flag = false;
  CoordSet* cs = nullptr;
  char ncolumn = 6;  // number of coordinates per line
  char nbyte = 12;   // width of one coordinate

  if (mode > 0) {
    ncolumn = 10;
    nbyte = 8;
  }
#define BUFSIZE 4194304
#define GETTING_LOW 10000

  else {
    if (I->CSTmpl) {
      cs = CoordSetCopy(I->CSTmpl);
    } else if (I->NCSet > 0) {
      cs = CoordSetCopy(I->CSet[0]);
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjMolLoadRSTFile: Missing topology" ENDFB(G);
      return (I);
    }
    CHECKOK(ok, cs);
    if (ok) {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " ObjMolLoadRSTFile: Loading from \"%s\".\n", fname ENDFB(G);
      p = buffer = FileGetContents(fname, nullptr);
      if (!buffer) ok = ErrMessage(G, __func__, "Unable to open file!");
    }
    if (ok) {
      p = nextline(p);
      if (mode == 0)  // skip NATOM,TIME
        p = nextline(p);
    }
    a = 0;
    b = 0;
    c = 0;
    f1 = 0.0;
    f2 = 0.0;
    while (ok && *p) {
      p = ncopy(cc, p, nbyte);
      if ((++b) == ncolumn) {
        b = 0;
        p = nextline(p);
      }
      f0 = f1;
      f1 = f2;
      if (sscanf(cc, "%f", &f2) == 1) {
        if ((++c) == 3) {
          c = 0;
          fp = cs->coordPtr(a);
          *(fp++) = f0;
          *(fp++) = f1;
          *(fp++) = f2;

          if ((++a) == I->NAtom) {
            a = 0;
            if (b) p = nextline(p);
            if (mode == 2)  // skip box
              p = nextline(p);
            b = 0;
            /* add new coord set */
            cs->invalidateRep(cRepAll, cRepInvRep);
            if (frame < 0) frame = I->NCSet;
            if (!I->NCSet) {
              zoom_flag = true;
            }

            VLACheck(I->CSet, CoordSet*, frame);
            CHECKOK(ok, I->CSet);
            if (ok) {
              if (I->NCSet <= frame) I->NCSet = frame + 1;
              delete I->CSet[frame];
              I->CSet[frame] = cs;
            }
            PRINTFB(G, FB_ObjectMolecule, FB_Details)
            " ObjectMolecule: read coordinates into state %d...\n",
                frame + 1 ENDFB(G);

            if (ok) cs = CoordSetCopy(cs);
            CHECKOK(ok, cs);

            if (mode == 0)  // restart file has only one frame
              break;

            frame += 1;
          }
        }
      } else {
        PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " ObjMolLoadRSTFile: atom/coordinate mismatch.\n" ENDFB(G);
        break;
      }
    }
    mfree(buffer);
  }
  delete cs;

  SceneChanged(G);
  SceneCountFrames(G);
  if (zoom_flag) {
    if (SettingGetGlobal_i(G, cSetting_auto_zoom)) {
      ExecutiveWindowZoom(G, I->Name, 0.0, -1, 0, 0,
                          quiet); /* auto zoom (all states) */
    }
  }
  return (I);
}

static const char* findflag(PyMOLGlobals* G, const char* p, const char* flag,
                            const char* format) {
  char cc[MAXLINELEN];
  char pat[MAXLINELEN] = "%FLAG ";
  int l;

  PRINTFD(G, FB_ObjectMolecule)
  " findflag: flag %s format %s\n", flag, format ENDFD;

  strcat(pat, flag);
  l = strlen(pat);
  while (*p) {
    p = ncopy(cc, p, l);
    if (WordMatch(G, cc, pat, true) < 0) {
      p = nextline(p);
      break;
    }
    p = nextline(p);
    if (!*p) {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjectMolecule-Error: Unrecognized file format (can't find \"%s\").\n",
          pat ENDFB(G);
    }
  }

  strcpy(pat, "%FORMAT(");
  strcat(pat, format);
  strcat(pat, ")");
  l = strlen(pat);
  while (*p) {
    p = ncopy(cc, p, l);
    if (WordMatch(G, cc, pat, true) < 0) {
      p = nextline(p);
      break;
    }
    p = nextline(p);
    if (!*p) {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjectMolecule-Error: Unrecognized file format (can't find \"%s\").\n",
          pat ENDFB(G);
    }
  }
  return (p);
}
static CoordSet* ObjectMoleculeTOPStr2CoordSet(PyMOLGlobals* G,
                                               const char* buffer,
                                               AtomInfoType** atInfoPtr) {
  const char* p;
  int nAtom;
  int a, b, c, bi, last_i, at_i, aa, rc;
  float* coord = nullptr;
  float* f;
  CoordSet* cset = nullptr;
  AtomInfoType *atInfo = nullptr, *ai;
  BondType *bond = nullptr, *bd;
  int nBond = 0;
  int auto_show = RepGetAutoShowMask(G);
  int amber7 = false;

  WordType title;
  ResName* resn;

  char cc[MAXLINELEN];
  int ok = true;
  int i0, i1, i2;

  /* trajectory parameters */

  int NTYPES, NBONH, MBONA, NTHETH, MTHETA;
  int NPHIH, MPHIA, NHPARM, NPARM, NNB, NRES;
  int NBONA, NTHETA, NPHIA, NUMBND, NUMANG, NPTRA;
  int NATYP, NPHB, IFPERT, NBPER, NGPER, NDPER;
  int MBPER, MGPER, MDPER, IFBOX = 0, NMXRS, IFCAP;
  int NEXTRA, IPOL = 0;
  int wid, col;
  float BETA;
  float BOX1, BOX2, BOX3;

  cset = CoordSetNew(G);

  p = buffer;
  nAtom = 0;
  if (atInfoPtr) atInfo = *atInfoPtr;
  if (!atInfo) ErrFatal(G, "TOPStr2CoordSet", "need atom information record!");
  /* failsafe for old version.. */

  ncopy(cc, p, 8);
  if (strcmp(cc, "%VERSION") == 0) {
    amber7 = true;
    PRINTFB(G, FB_ObjectMolecule, FB_Details)
    " ObjectMolecule: Attempting to read Amber7 topology file.\n" ENDFB(G);
  } else {
    PRINTFB(G, FB_ObjectMolecule, FB_Details)
    " ObjectMolecule: Assuming this is an Amber6 topology file.\n" ENDFB(G);
  }

  /* read title */
  if (amber7) {
    p = findflag(G, buffer, "TITLE", "20a4");
  }

  p = ncopy(cc, p, 20);
  title[0] = 0;
  sscanf(cc, "%s", title);
  p = nextline(p);

  if (amber7) {
    p = findflag(G, buffer, "POINTERS", "10I8");

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &nAtom) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NTYPES) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NBONH) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MBONA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NTHETH) == 1);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MTHETA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPHIH) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MPHIA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NHPARM) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPARM) == 1);

    p = nextline(p);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NNB) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NRES) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NBONA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NTHETA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPHIA) == 1);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NUMBND) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NUMANG) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPTRA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NATYP) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPHB) == 1);

    p = nextline(p);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &IFPERT) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NBPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NGPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NDPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MBPER) == 1);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MGPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MDPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &IFBOX) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NMXRS) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &IFCAP) == 1);

    p = nextline(p);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NEXTRA) == 1);

  } else {
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &nAtom) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NTYPES) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NBONH) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MBONA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NTHETH) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MTHETA) == 1);

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPHIH) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MPHIA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NHPARM) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPARM) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NNB) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NRES) == 1);

    p = nextline(p);

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NBONA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NTHETA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPHIA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NUMBND) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NUMANG) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPTRA) == 1);

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NATYP) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPHB) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &IFPERT) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NBPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NGPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NDPER) == 1);

    p = nextline(p);

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MBPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MGPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MDPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &IFBOX) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NMXRS) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &IFCAP) == 1);

    p = ncopy(cc, p, 6);
    if (sscanf(cc, "%d", &NEXTRA) != 1) NEXTRA = 0;
  }

  if (!ok) {
    ErrMessage(G, "TOPStrToCoordSet", "Error reading counts lines");
    ok_raise(1);
  } else {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    " TOPStr2CoordSet: read counts line nAtom %d NBONA %d NBONH %d\n", nAtom,
        NBONA, NBONH ENDFB(G);
  }

  switch (IFBOX) {
    case 2:
      cset->PeriodicBoxType = CoordSet::Octahedral;
      PRINTFB(G, FB_ObjectMolecule, FB_Details)
      " TOPStrToCoordSet: Warning: can't currently image a truncated "
      "octahedron...\n" ENDFB(G);
      break;
    case 1:
      cset->PeriodicBoxType = CoordSet::Orthogonal;
      break;
    case 0:
    default:
      cset->PeriodicBoxType = CoordSet::NoPeriodicity;
      break;
  }

  p = nextline(p);

  if (ok) {
    VLACheck(atInfo, AtomInfoType, nAtom);

    if (amber7) {
      p = findflag(G, buffer, "ATOM_NAME", "20a4");
    }
    /* read atoms */

    b = 0;
    for (a = 0; a < nAtom; a++) {
      p = ntrim(cc, p, 4);
      atInfo[a].name = LexIdx(G, cc);
      if ((++b) == 20) {
        b = 0;
        p = nextline(p);
      }
    }

    if (b) p = nextline(p);

    if (!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading atom names");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read atom names.\n" ENDFB(G);
    }

    /* read charges */

    if (amber7) {
      p = findflag(G, buffer, "CHARGE", "5E16.8");
    }

    b = 0;
    for (a = 0; a < nAtom; a++) {
      p = ncopy(cc, p, 16);
      ai = atInfo + a;
      if (!sscanf(cc, "%f", &ai->partialCharge))
        ok = false;
      else {
        ai->partialCharge /= 18.2223F; /* convert to electron charge */
      }
      if ((++b) == 5) {
        b = 0;
        p = nextline(p);
      }
    }

    if (!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading charges");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read charges.\n" ENDFB(G);
    }
    if (b) p = nextline(p);

    if (!amber7) {
      /* skip masses */

      p = skip_fortran(nAtom, 5, p);
    }

    /* read LJ atom types */

    if (amber7) {
      p = findflag(G, buffer, "ATOM_TYPE_INDEX", "10I8");
      col = 10;
      wid = 8;
    } else {
      col = 12;
      wid = 6;
    }

    b = 0;
    for (a = 0; a < nAtom; a++) {
      p = ncopy(cc, p, wid);
      ai = atInfo + a;
      if (!sscanf(cc, "%d", &ai->customType)) ok = false;
      if ((++b) == col) {
        b = 0;
        p = nextline(p);
      }
    }
    if (b) p = nextline(p);

    if (!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error LJ atom types");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read LJ atom types.\n" ENDFB(G);
    }

    if (!amber7) {
      /* skip excluded atom counts */

      p = skip_fortran(nAtom, 12, p);

      /* skip NB param arrays */

      p = skip_fortran(NTYPES * NTYPES, 12, p);
    }

    /* read residue labels */

    if (amber7) {
      p = findflag(G, buffer, "RESIDUE_LABEL", "20a4");
    }

    resn = pymol::malloc<ResName>(NRES);

    b = 0;
    for (a = 0; a < NRES; a++) {
      p = ncopy(cc, p, 4);
      if (!sscanf(cc, "%s", resn[a])) resn[a][0] = 0;
      if ((++b) == 20) {
        b = 0;
        p = nextline(p);
      }
    }
    if (b) p = nextline(p);

    if (!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading residue labels");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read residue labels.\n" ENDFB(G);
    }

    /* read residue assignments */

    if (amber7) {
      p = findflag(G, buffer, "RESIDUE_POINTER", "10I8");
      col = 10;
      wid = 8;
    } else {
      col = 12;
      wid = 6;
    }

    b = 0;
    last_i = 0;
    rc = 0;
    for (a = 0; a < NRES; a++) {
      p = ncopy(cc, p, wid);
      if (sscanf(cc, "%d", &at_i)) {
        if (last_i)
          for (aa = (last_i - 1); aa < (at_i - 1); aa++) {
            ai = atInfo + aa;
            ai->resn = LexIdx(G, resn[a - 1]);
            ai->resv = rc;
          }
        rc++;
        last_i = at_i;
      }
      if ((++b) == col) {
        b = 0;
        p = nextline(p);
      }
    }
    if (b) p = nextline(p);
    if (last_i)
      for (aa = (last_i - 1); aa < nAtom; aa++) {
        ai = atInfo + aa;
        ai->resn = LexIdx(G, resn[NRES - 1]);
        ai->resv = rc;
      }
    rc++;

    if (!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading residues");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read residues.\n" ENDFB(G);
    }

    FreeP(resn);

    if (!amber7) {
      /* skip bond force constants */

      p = skip_fortran(NUMBND, 5, p);

      /* skip bond lengths */

      p = skip_fortran(NUMBND, 5, p);

      /* skip angle force constant */

      p = skip_fortran(NUMANG, 5, p);

      /* skip angle eq */

      p = skip_fortran(NUMANG, 5, p);

      /* skip dihedral force constant */

      p = skip_fortran(NPTRA, 5, p);

      /* skip dihedral periodicity */

      p = skip_fortran(NPTRA, 5, p);

      /* skip dihedral phases */

      p = skip_fortran(NPTRA, 5, p);

      /* skip SOLTYs */

      p = skip_fortran(NATYP, 5, p);

      /* skip LJ terms r12 */

      p = skip_fortran((NTYPES * (NTYPES + 1)) / 2, 5, p);

      /* skip LJ terms r6 */

      p = skip_fortran((NTYPES * (NTYPES + 1)) / 2, 5, p);
    }

    /* read bonds */

    if (amber7) {
      p = findflag(G, buffer, "BONDS_INC_HYDROGEN", "10I8");
      col = 10;
      wid = 8;
    } else {
      col = 12;
      wid = 6;
    }

    nBond = NBONH + NBONA;

    bond = VLACalloc(BondType, nBond);

    bi = 0;

    b = 0;
    c = 0;
    i0 = 0;
    i1 = 0;
    for (a = 0; a < 3 * NBONH; a++) {
      p = ncopy(cc, p, wid);
      i2 = i1;
      i1 = i0;
      if (!sscanf(cc, "%d", &i0)) ok = false;
      if ((++c) == 3) {
        c = 0;
        bd = bond + bi;
        bd->index[0] = (abs(i2) / 3);
        bd->index[1] = (abs(i1) / 3);
        bd->order = 1;
        bi++;
      }
      if ((++b) == col) {
        b = 0;
        p = nextline(p);
      }
    }
    if (b) p = nextline(p);

    if (!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error hydrogen containing bonds");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read %d hydrogen containing bonds.\n", NBONH ENDFB(G);
    }

    if (amber7) {
      p = findflag(G, buffer, "BONDS_WITHOUT_HYDROGEN", "10I8");
      col = 10;
      wid = 8;
    } else {
      col = 12;
      wid = 6;
    }

    b = 0;
    c = 0;
    for (a = 0; a < 3 * NBONA; a++) {
      p = ncopy(cc, p, wid);
      i2 = i1;
      i1 = i0;
      if (!sscanf(cc, "%d", &i0)) ok = false;
      if ((++c) == 3) {
        c = 0;
        bd = bond + bi;
        bd->index[0] = (abs(i2) / 3);
        bd->index[1] = (abs(i1) / 3);
        bd->order = 1;  // PYMOL-2707
        bi++;
      }
      if ((++b) == col) {
        b = 0;
        p = nextline(p);
      }
    }
    if (b) p = nextline(p);

    if (!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error hydrogen free bonds");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read %d hydrogen free bonds.\n", NBONA ENDFB(G);
    }

    if (!amber7) {
      /* skip hydrogen angles */

      p = skip_fortran(4 * NTHETH, 12, p);

      /* skip non-hydrogen angles */

      p = skip_fortran(4 * NTHETA, 12, p);

      /* skip hydrogen dihedrals */

      p = skip_fortran(5 * NPHIH, 12, p);

      /* skip non hydrogen dihedrals */

      p = skip_fortran(5 * NPHIA, 12, p);

      /* skip nonbonded exclusions */

      p = skip_fortran(NNB, 12, p);

      /* skip hydrogen bonds ASOL */

      p = skip_fortran(NPHB, 5, p);

      /* skip hydrogen bonds BSOL */

      p = skip_fortran(NPHB, 5, p);

      /* skip HBCUT */

      p = skip_fortran(NPHB, 5, p);
    }
    /* read AMBER atom types */

    if (amber7) {
      p = findflag(G, buffer, "AMBER_ATOM_TYPE", "20a4");
    }

    b = 0;
    for (a = 0; a < nAtom; a++) {
      OrthoLineType temp;
      p = ncopy(cc, p, 4);
      ai = atInfo + a;
      if (sscanf(cc, "%s", temp) != 1)
        ok = false;
      else {
        ai->textType = LexIdx(G, temp);
      }
      if ((++b) == 20) {
        b = 0;
        p = nextline(p);
      }
    }
    if (b) p = nextline(p);

    if (!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading atom types");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read atom types.\n" ENDFB(G);
    }

    if (!amber7) {
      /* skip TREE classification */

      p = skip_fortran(nAtom, 20, p);

      /* skip tree joining information */

      p = skip_fortran(nAtom, 12, p);

      /* skip last atom rotated blah blah blah */

      p = skip_fortran(nAtom, 12, p);
    }

    if (IFBOX > 0) {
      int IPTRES, NSPM = 0, NSPSOL;

      if (amber7) {
        p = findflag(G, buffer, "SOLVENT_POINTERS", "3I8");
        wid = 8;
      } else {
        wid = 6;
      }
      p = ncopy(cc, p, wid);
      ok = ok && (sscanf(cc, "%d", &IPTRES) == 1);
      p = ncopy(cc, p, wid);
      ok = ok && (sscanf(cc, "%d", &NSPM) == 1);
      p = ncopy(cc, p, wid);
      ok = ok && (sscanf(cc, "%d", &NSPSOL) == 1);

      ok_assert(1, ok);

      p = nextline(p);

      if (amber7) {
        p = findflag(G, buffer, "ATOMS_PER_MOLECULE", "10I8");
        col = 10;
      } else {
        col = 12;
      }

      /* skip num atoms per box */

      p = skip_fortran(NSPM, col, p);

      if (amber7) {
        p = findflag(G, buffer, "BOX_DIMENSIONS", "5E16.8");
      }
      wid = 16;

      p = ncopy(cc, p, 16);
      ok = ok && (sscanf(cc, "%f", &BETA) == 1);
      p = ncopy(cc, p, 16);
      ok = ok && (sscanf(cc, "%f", &BOX1) == 1);
      p = ncopy(cc, p, 16);
      ok = ok && (sscanf(cc, "%f", &BOX2) == 1);
      p = ncopy(cc, p, 16);
      ok = ok && (sscanf(cc, "%f", &BOX3) == 1);

      if (ok) {
        float angle[3] = {90.f, BETA, 90.f};

        if ((BETA > 109.47) && (BETA < 109.48)) {
          cset->PeriodicBoxType = CoordSet::Octahedral;
          angle[0] = 109.47122;
          angle[1] = 109.47122;
          angle[2] = 109.47122;
        } else if (BETA == 60.0) {
          angle[0] = 60.0; /* rhombic dodecahedron (from ptraj.c) */
          angle[1] = 90.0;
          angle[2] = 60.0;
        }

        cset->Symmetry = std::make_unique<CSymmetry>(G);
        cset->Symmetry->Crystal.setDims(BOX1, BOX2, BOX3);
        cset->Symmetry->Crystal.setAngles(angle);
      }
      /* skip periodic box */

      p = nextline(p);
    }

    if (!amber7) {
      if (IFCAP > 0) {
        p = nextline(p);
        p = nextline(p);
        p = nextline(p);
      }

      if (IFPERT > 0) {
        /* skip perturbed bond atoms */

        p = skip_fortran(2 * NBPER, 12, p);

        /* skip perturbed bond atom pointers */

        p = skip_fortran(2 * NBPER, 12, p);

        /* skip perturbed angles */

        p = skip_fortran(3 * NGPER, 12, p);

        /* skip perturbed angle pointers */

        p = skip_fortran(2 * NGPER, 12, p);

        /* skip perturbed dihedrals */

        p = skip_fortran(4 * NDPER, 12, p);

        /* skip perturbed dihedral pointers */

        p = skip_fortran(2 * NDPER, 12, p);

        /* skip residue names */

        p = skip_fortran(NRES, 20, p);

        /* skip atom names */

        p = skip_fortran(nAtom, 20, p);

        /* skip atom symbols */

        p = skip_fortran(nAtom, 20, p);

        /* skip unused field */

        p = skip_fortran(nAtom, 5, p);

        /* skip perturbed flags */

        p = skip_fortran(nAtom, 12, p);

        /* skip LJ atom flags */

        p = skip_fortran(nAtom, 12, p);

        /* skip perturbed charges */

        p = skip_fortran(nAtom, 5, p);
      }

      if (IPOL > 0) {
        /* skip atomic polarizabilities */

        p = skip_fortran(nAtom, 5, p);
      }

      if ((IPOL > 0) && (IFPERT > 0)) {
        /* skip atomic polarizabilities */

        p = skip_fortran(nAtom, 5, p);
      }
    }
    /* for future reference

       %FLAG LES_NTYP
       %FORMAT(10I8)
       %FLAG LES_TYPE
       %FORMAT(10I8)
       %FLAG LES_FAC
       %FORMAT(5E16.8)
       %FLAG LES_CNUM
       %FORMAT(10I8)
       %FLAG LES_ID
       %FORMAT(10I8)

       Here is the additional information for LES topology formats:
       First, if NPARM ==1, LES entries are in topology (NPARM is the 10th
       entry in the initial list of control parameters); otherwise the standard
       format applies.
       So, with NPARM=1, you just need to read a few more things at the very
       end of topology file:
       LES_NTYP (format: I6) ... one number, number of LES types
       and four arrays:
       LES_TYPE (12I6) ... NATOM integer entries
       LES_FAC (E16.8) ... LES_NTYPxLES_NTYP float entries
       LES_CNUM (12I6) ... NATOM integer entries
       LES_ID (12I6)   ... NATOM integer entries

       and that's it. Your parser must have skipped this information because it
       was at the end of the file. Maybe that's good enough.

     */

    coord = VLAlloc(float, 3 * nAtom);

    f = coord;
    for (a = 0; a < nAtom; a++) {
      *(f++) = 0.0;
      *(f++) = 0.0;
      *(f++) = 0.0;
      ai = atInfo + a;
      ai->id = a + 1; /* assign 1-based identifiers */
      ai->rank = a;
      AtomInfoAssignParameters(G, ai);
      AtomInfoAssignColors(G, ai);
      ai->visRep = auto_show;
    }
  }
  if (ok) {
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->TmpBond = pymol::vla_take_ownership(bond);
    cset->NTmpBond = nBond;
  } else {
  ok_except1:
    delete cset;
    cset = nullptr;
    ErrMessage(G, __func__, "failed");
  }
  if (atInfoPtr) *atInfoPtr = atInfo;

  return (cset);
}
static ObjectMolecule* ObjectMoleculeReadTOPStr(PyMOLGlobals* G,
                                                ObjectMolecule* I, char* TOPStr,
                                                int frame, int discrete) {
  CoordSet* cset = nullptr;
  pymol::vla<AtomInfoType> atInfo(1);
  int ok = true;
  int isNew = true;
  unsigned int nAtom = 0;

  if (!I)
    isNew = true;
  else
    isNew = false;

  if (ok) {
    if (isNew) {
      I = (ObjectMolecule*)new ObjectMolecule(G, discrete);
      CHECKOK(ok, I);
      if (ok) std::swap(atInfo, I->AtomInfo);
    }
    if (ok && isNew) {
      I->Color = AtomInfoUpdateAutoColor(G);
    }

    if (ok) cset = ObjectMoleculeTOPStr2CoordSet(G, TOPStr, &atInfo);
    CHECKOK(ok, cset);
  }

  /* include coordinate set */
  if (ok) {
    nAtom = cset->NIndex;

    if (I->DiscreteFlag && atInfo) {
      unsigned int a;
      int fp1 = frame + 1;
      AtomInfoType* ai = atInfo.data();
      for (a = 0; a < nAtom; a++) {
        (ai++)->discrete_state = fp1;
      }
    }

    cset->Obj = I;
    cset->enumIndices();
    cset->invalidateRep(cRepAll, cRepInvRep);
    if (isNew) {
      std::swap(I->AtomInfo, atInfo);
    } else if (ok) {
      ok &= ObjectMoleculeMerge(I, std::move(atInfo), cset, false, cAIC_AllMask,
                                true);
    }
    if (isNew) I->NAtom = nAtom;
    /*
       if(frame<0) frame=I->NCSet;
       VLACheck(I->CSet,CoordSet*,frame);
       if(I->NCSet<=frame) I->NCSet=frame+1;
       if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
       I->CSet[frame] = cset;
     */

    if (ok && isNew) ok = ObjectMoleculeConnect(I, cset, false);
    if (cset->Symmetry && (!I->Symmetry)) {
      I->Symmetry.reset(new CSymmetry(*cset->Symmetry));
      CHECKOK(ok, I->Symmetry);
    }

    delete I->CSTmpl;
    I->CSTmpl = cset; /* save template coordinate set */

    SceneCountFrames(G);
    if (ok) ok &= ObjectMoleculeExtendIndices(I, -1);
    if (ok) ok &= ObjectMoleculeSort(I);
    if (ok) {
      ObjectMoleculeUpdateIDNumbers(I);
      ObjectMoleculeUpdateNonbonded(I);
    }
  }
  if (!ok) {
    DeleteP(I)
  }
  return (I);
}

ObjectMolecule* ObjectMoleculeLoadTOPFile(PyMOLGlobals* G, ObjectMolecule* obj,
                                          const char* fname, int frame,
                                          int discrete) {
  ObjectMolecule* I = nullptr;
  char* buffer;

  buffer = FileGetContents(fname, nullptr);

  if (!buffer)
    ErrMessage(G, __func__, "Unable to open file!");
  else {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    " %s: Loading from %s.\n", __func__, fname ENDFB(G);

    I = ObjectMoleculeReadTOPStr(G, obj, buffer, frame, discrete);
    mfree(buffer);
  }

  return (I);
}
