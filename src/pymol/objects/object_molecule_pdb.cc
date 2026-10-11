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

// Reading PDB and PQR files.

#include "pymol/base/parse.h"
#include "pymol/base/err.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/base/feedback.h"
#include "pymol/base/util2.h"
#include "pymol/executive/executive.h"
#include "pymol/objects/object_molecule_internal.h"

#define ntrim ParseNTrim
#define nextline ParseNextLine
#define ncopy ParseNCopy
#define nskip ParseNSkip

#define cResvMask 0x7FFF

#define strstartswith p_strstartswith

static int strstartswithword(const char* s, const char* word) {
  while (*word)
    if (*s++ != *word++) return false;
  switch (*s) {
    case ' ':
    case '\t':
    case '\r':
    case '\n':
    case '\0':
      return true;
  }
  return false;
}

static void AddOrthoOutputIfMatchesTags(PyMOLGlobals* G, int n_tags, int nAtom,
                                        const char* const* tag_start,
                                        const char* p, char* cc, int quiet) {
  if (n_tags && !quiet && !(nAtom > 0 && strstartswith(p, "HEADER"))) {
    // HEADER is the mark for a new object, this is a new object, and
    // gets processed on the next pass, when nAtom=0
    int tc = 0;
    for (; tc < n_tags; tc++) {
      if (!strstartswithword(p, tag_start[tc])) continue;
      ParseNTrimRight(cc, p, MAXLINELEN - 1);
      OrthoAddOutput(G, cc);
      OrthoNewLine(G, nullptr, true);
      break;
    }
  }
}

/**
 * Based on PDB nomenclature (resn, name), do:
 *
 * 1) If `ai1` or `ai2` is a known charged PDB atom, assign `formalCharge`
 *    and set `chemFlag` to false
 *
 * 2) If `ai1` and `ai2` are connected by a double bond, set (*bond_order) = 2
 */
void assign_pdb_known_residue(PyMOLGlobals* G, AtomInfoType* ai1,
                              AtomInfoType* ai2, int* bond_order) {
  int order = *(bond_order);
  auto const name1 = pymol::zstring_view(LexStr(G, ai1->name));
  auto const name2 = pymol::zstring_view(LexStr(G, ai2->name));
  auto const resn1 = pymol::zstring_view(LexStr(G, ai1->resn));

  /* nasty high-speed hack to get bond valences and formal charges
     for standard residues */
  if (((name1 == "C" && name2 == "O") || (name2 == "C" && name1 == "O")) &&
      AtomInfoKnownProteinResName(resn1.c_str())) {
    order = 2;
  } else if (name2 == "C" && name1 == "OXT") {
    ai1->formalCharge = -1;
    ai1->chemFlag = false;
  } else if (name1 == "C" && name2 == "OXT") {
    ai2->formalCharge = -1;
    ai2->chemFlag = false;
  } else {
    switch (resn1[0]) {
      case 'A':
        switch (resn1[1]) {
          case 'R':
            switch (resn1[2]) {
              case 'G': /* ARG... */
                switch (resn1[3]) {
                  case 0:
                  case 'P': /*  ARG, ARGP */
                    if (name1 == "NH1") {
                      ai1->formalCharge = 1;
                      ai1->chemFlag = false;
                    } else if (name2 == "NH1") {
                      ai2->formalCharge = 1;
                      ai2->chemFlag = false;
                    }
                    // PYMOL-5019: If we're forcing NH1 to be formal charge 1,
                    // the other atom must have formal charge 0.
                    if (name1 == "NH2") {
                      ai1->formalCharge = 0;
                      ai1->chemFlag = false;
                    } else if (name2 == "NH2") {
                      ai2->formalCharge = 0;
                      ai2->chemFlag = false;
                    }
                    break;
                }
                if ((name1 == "CZ" && name2 == "NH1") ||
                    (name2 == "CZ" && name1 == "NH1"))
                  order = 2;
                // PYMOL-5019: If we're forcing NH1 to be double bounded, the
                // other bond must display as single.
                if ((name1 == "CZ" && name2 == "NH2") ||
                    (name2 == "CZ" && name1 == "NH2"))
                  order = 1;
                break;
            }
            break;
          case 'S':
            switch (resn1[2]) {
              case 'P': /* ASP... */
                switch (resn1[3]) {
                  case 0:
                  case 'M': /* ASP, ASPM minus assumption */
                    if (name1 == "OD2") {
                      ai1->formalCharge = -1;
                      ai1->chemFlag = false;
                    } else if (name2 == "OD2") {
                      ai2->formalCharge = -1;
                      ai2->chemFlag = false;
                    }
                    break;
                }
                if ((name1 == "CG" && name2 == "OD1") ||
                    (name2 == "CG" && name1 == "OD1"))
                  order = 2;
                break;
              case 'N': /* ASN  */
                if ((name1 == "CG" && name2 == "OD1") ||
                    (name2 == "CG" && name1 == "OD1"))
                  order = 2;
                break;
            }
            break;
          case 0:
            if ((name1 == "O2P" || name1 == "OP2")) {
              ai1->formalCharge = -1;
              ai1->chemFlag = false;
            } else if ((name2 == "O2P" || name2 == "OP2")) {
              ai2->formalCharge = -1;
              ai2->chemFlag = false;
            }
            if ((name1 == "C8" && name2 == "N7") ||
                (name2 == "C8" && name1 == "N7"))
              order = 2;
            else if ((name1 == "C4" && name2 == "C5") ||
                     (name2 == "C4" && name1 == "C5"))
              order = 2;

            else if ((name1 == "C6" && name2 == "N1") ||
                     (name2 == "C6" && name1 == "N1"))
              order = 2;
            else if ((name1 == "C2" && name2 == "N3") ||
                     (name2 == "C2" && name1 == "N3"))
              order = 2;
            else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                     (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
              order = 2;
            break;
        }
        break;
      case 'C':
        if (resn1[1] == 0) {
          if ((name1 == "O2P" || name1 == "OP2")) {
            ai1->formalCharge = -1;
            ai1->chemFlag = false;
          } else if ((name2 == "O2P" || name2 == "OP2")) {
            ai2->formalCharge = -1;
            ai2->chemFlag = false;
          }
          if ((name1 == "C2" && name2 == "O2") ||
              (name2 == "C2" && name1 == "O2"))
            order = 2;
          else if ((name1 == "C4" && name2 == "N3") ||
                   (name2 == "C4" && name1 == "N3"))
            order = 2;

          else if ((name1 == "C5" && name2 == "C6") ||
                   (name2 == "C5" && name1 == "C6"))
            order = 2;
          else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                   (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
            order = 2;
        }
        break;
      case 'D': /* deoxy nucleic acids */
        switch (resn1[1]) {
          case 'A':
            if (resn1[2] == 0) {
              if ((name1 == "O2P" || name1 == "OP2")) {
                ai1->formalCharge = -1;
                ai1->chemFlag = false;
              } else if ((name2 == "O2P" || name2 == "OP2")) {
                ai2->formalCharge = -1;
                ai2->chemFlag = false;
              }
              if ((name1 == "C8" && name2 == "N7") ||
                  (name2 == "C8" && name1 == "N7"))
                order = 2;
              else if ((name1 == "C4" && name2 == "C5") ||
                       (name2 == "C4" && name1 == "C5"))
                order = 2;

              else if ((name1 == "C6" && name2 == "N1") ||
                       (name2 == "C6" && name1 == "N1"))
                order = 2;
              else if ((name1 == "C2" && name2 == "N3") ||
                       (name2 == "C2" && name1 == "N3"))
                order = 2;
              else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                       (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
                order = 2;
            }
            break;
          case 'C':
            if (resn1[2] == 0) {
              if ((name1 == "O2P" || name1 == "OP2")) {
                ai1->formalCharge = -1;
                ai1->chemFlag = false;
              } else if ((name2 == "O2P" || name2 == "OP2")) {
                ai2->formalCharge = -1;
                ai2->chemFlag = false;
              }
              if ((name1 == "C2" && name2 == "O2") ||
                  (name2 == "C2" && name1 == "O2"))
                order = 2;
              else if ((name1 == "C4" && name2 == "N3") ||
                       (name2 == "C4" && name1 == "N3"))
                order = 2;

              else if ((name1 == "C5" && name2 == "C6") ||
                       (name2 == "C5" && name1 == "C6"))
                order = 2;
              else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                       (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
                order = 2;
            }
            break;
          case 'T':
            if (resn1[2] == 0) {
              if ((name1 == "O2P" || name1 == "OP2"))
                ai1->formalCharge = -1;
              else if ((name2 == "O2P" || name2 == "OP2"))
                ai2->formalCharge = -1;

              if ((name1 == "C2" && name2 == "O2") ||
                  (name2 == "C2" && name1 == "O2"))
                order = 2;
              else if ((name1 == "C4" && name2 == "O4") ||
                       (name2 == "C4" && name1 == "O4"))
                order = 2;

              else if ((name1 == "C5" && name2 == "C6") ||
                       (name2 == "C5" && name1 == "C6"))
                order = 2;
              else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                       (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
                order = 2;
            }
            break;
          case 'G':
            if (resn1[2] == 0) {
              if ((name1 == "O2P" || name1 == "OP2")) {
                ai1->formalCharge = -1;
                ai1->chemFlag = false;
              } else if ((name2 == "O2P" || name2 == "OP2")) {
                ai2->formalCharge = -1;
                ai2->chemFlag = false;
              }
              if ((name1 == "C6" && name2 == "O6") ||
                  (name2 == "C6" && name1 == "O6"))
                order = 2;
              else if ((name1 == "C2" && name2 == "N3") ||
                       (name2 == "C2" && name1 == "N3"))
                order = 2;
              else if ((name1 == "C8" && name2 == "N7") ||
                       (name2 == "C8" && name1 == "N7"))
                order = 2;
              else if ((name1 == "C4" && name2 == "C5") ||
                       (name2 == "C4" && name1 == "C5"))
                order = 2;
              else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                       (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
                order = 2;
            }
            break;
          case 'U':
            if (resn1[2] == 0) {
              if ((name1 == "O2P" || name1 == "OP2")) {
                ai1->formalCharge = -1;
                ai1->chemFlag = false;
              } else if ((name2 == "O2P" || name2 == "OP2")) {
                ai2->formalCharge = -1;
                ai2->chemFlag = false;
              }

              if ((name1 == "C2" && name2 == "O2") ||
                  (name2 == "C2" && name1 == "O2"))
                order = 2;
              else if ((name1 == "C4" && name2 == "O4") ||
                       (name2 == "C4" && name1 == "O4"))
                order = 2;

              else if ((name1 == "C5" && name2 == "C6") ||
                       (name2 == "C5" && name1 == "C6"))
                order = 2;
              else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                       (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
                order = 2;
            }
            break;
        }
        break;
      case 'G':
        switch (resn1[1]) {
          case 'L':
            switch (resn1[2]) {
              case 'U': /* GLU missing GLUN, GLUH, GLH handling */
                switch (resn1[3]) {
                  case 0:
                  case 'M': /* minus */
                    if (name1 == "OE2") {
                      ai1->formalCharge = -1;
                      ai1->chemFlag = false;
                    } else if (name2 == "OE2") {
                      ai2->formalCharge = -1;
                      ai2->chemFlag = false;
                    }
                    break;
                }
                if ((name1 == "CD" && name2 == "OE1") ||
                    (name2 == "CD" && name1 == "OE1"))
                  order = 2;
                break;
              case 'N': /* GLN or GLU */
                if ((name1 == "CD" && name2 == "OE1") ||
                    (name2 == "CD" && name1 == "OE1"))
                  order = 2;
                break;
            }
            break;
          case 0:
            if ((name1 == "O2P" || name1 == "OP2")) {
              ai1->formalCharge = -1;
              ai1->chemFlag = false;
            } else if ((name2 == "O2P" || name2 == "OP2")) {
              ai2->formalCharge = -1;
              ai2->chemFlag = false;
            }

            if ((name1 == "C6" && name2 == "O6") ||
                (name2 == "C6" && name1 == "O6"))
              order = 2;
            else if ((name1 == "C2" && name2 == "N3") ||
                     (name2 == "C2" && name1 == "N3"))
              order = 2;
            else if ((name1 == "C8" && name2 == "N7") ||
                     (name2 == "C8" && name1 == "N7"))
              order = 2;
            else if ((name1 == "C4" && name2 == "C5") ||
                     (name2 == "C4" && name1 == "C5"))
              order = 2;
            else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                     (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
              order = 2;
            break;
        }
        break;
      case 'H':
        switch (resn1[1]) {
          case 'I':
            switch (resn1[2]) {
              case 'P':
                if (name1 == "ND1") {
                  ai1->formalCharge = 1;
                  ai1->chemFlag = false;
                } else if (name2 == "ND1") {
                  ai2->formalCharge = 1;
                  ai2->chemFlag = false;
                }
                if ((name1 == "CG" && name2 == "CD2") ||
                    (name2 == "CG" && name1 == "CD2"))
                  order = 2;
                else if ((name1 == "CE1" && name2 == "ND1") ||
                         (name2 == "CE1" && name1 == "ND1"))
                  order = 2;
                break;
              case 'S':
                switch (resn1[3]) {
                  case 'A': /* HISA Gromacs */
                  case 'D': /* HISD Quanta */
                    if ((name1 == "CG" && name2 == "CD2") ||
                        (name2 == "CG" && name1 == "CD2"))
                      order = 2;
                    else if ((name1 == "CE1" && name2 == "NE2") ||
                             (name2 == "CE1" && name1 == "NE2"))
                      order = 2;
                    break;
                  case 0:   /* plain HIS */
                  case 'B': /* HISB Gromacs */
                  case 'E': /* HISE Quanta */
                    if ((name1 == "CG" && name2 == "CD2") ||
                        (name2 == "CG" && name1 == "CD2"))
                      order = 2;
                    else if ((name1 == "CE1" && name2 == "ND1") ||
                             (name2 == "CE1" && name1 == "ND1"))
                      order = 2;
                    break;
                  case 'H': /* HISH Gromacs */
                  case 'P': /* HISP Quanta */
                    if (name1 == "ND1") {
                      ai1->formalCharge = 1;
                      ai1->chemFlag = false;
                    } else if (name2 == "ND1") {
                      ai2->formalCharge = 1;
                      ai2->chemFlag = false;
                    }
                    if ((name1 == "CG" && name2 == "CD2") ||
                        (name2 == "CG" && name1 == "CD2"))
                      order = 2;
                    else if ((name1 == "CE1" && name2 == "ND1") ||
                             (name2 == "CE1" && name1 == "ND1"))
                      order = 2;
                    break;
                }
                break;
              case 'E': /* HIE */
                if ((name1 == "CG" && name2 == "CD2") ||
                    (name2 == "CG" && name1 == "CD2"))
                  order = 2;
                else if ((name1 == "CE1" && name2 == "ND1") ||
                         (name2 == "CE1" && name1 == "ND1"))
                  order = 2;
                break;
              case 'D': /* HID */
                if ((name1 == "CG" && name2 == "CD2") ||
                    (name2 == "CG" && name1 == "CD2"))
                  order = 2;
                else if ((name1 == "CE1" && name2 == "NE2") ||
                         (name2 == "CE1" && name1 == "NE2"))
                  order = 2;
                break;
            }
            break;
        }
        break;
      case 'I':
        if (resn1[1] == 0) {
          if ((name1 == "O2P" || name1 == "OP2")) {
            ai1->formalCharge = -1;
            ai1->chemFlag = false;
          } else if ((name2 == "O2P" || name2 == "OP2")) {
            ai2->formalCharge = -1;
            ai2->chemFlag = false;
          }
          if ((name1 == "C8" && name2 == "N7") ||
              (name2 == "C8" && name1 == "N7"))
            order = 2;
          else if ((name1 == "C4" && name2 == "C5") ||
                   (name2 == "C4" && name1 == "C5"))
            order = 2;

          else if ((name1 == "C6" && name2 == "N1") ||
                   (name2 == "C6" && name1 == "N1"))
            order = 2;
          else if ((name1 == "C2" && name2 == "N3") ||
                   (name2 == "C2" && name1 == "N3"))
            order = 2;
          else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                   (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
            order = 2;
        }
        break;
      case 'P':
        switch (resn1[1]) {
          case 'H': /* PHE */
            if (resn1[2] == 'E') {
              if ((name1 == "CG" && name2 == "CD1") ||
                  (name2 == "CG" && name1 == "CD1"))
                order = 2;
              else if ((name1 == "CZ" && name2 == "CE1") ||
                       (name2 == "CZ" && name1 == "CE1"))
                order = 2;

              else if ((name1 == "CE2" && name2 == "CD2") ||
                       (name2 == "CE2" && name1 == "CD2"))
                order = 2;
              break;
            }
        }
        break;
      case 'L':
        switch (resn1[1]) {
          case 'Y':
            switch (resn1[2]) {
              case 'S': /* LYS. */
                switch (resn1[3]) {
                  case 0:
                  case 'P': /* LYS, LYSP */
                    if (name1 == "NZ") {
                      ai1->formalCharge = 1;
                      ai1->chemFlag = false;
                    } else if (name2 == "NZ") {
                      ai2->formalCharge = 1;
                      ai2->chemFlag = false;
                    }
                    break;
                }
                break;
            }
            break;
        }
        break;
      case 'T':
        switch (resn1[1]) {
          case 'Y': /* TYR */
            if (resn1[2] == 'R') {
              if ((name1 == "CG" && name2 == "CD1") ||
                  (name2 == "CG" && name1 == "CD1"))
                order = 2;
              else if ((name1 == "CZ" && name2 == "CE1") ||
                       (name2 == "CZ" && name1 == "CE1"))
                order = 2;

              else if ((name1 == "CE2" && name2 == "CD2") ||
                       (name2 == "CE2" && name1 == "CD2"))
                order = 2;
              break;
            }
            break;
          case 'R':
            if (resn1[2] == 'P') {
              if ((name1 == "CG" && name2 == "CD1") ||
                  (name2 == "CG" && name1 == "CD1"))
                order = 2;
              else if ((name1 == "CZ3" && name2 == "CE3") ||
                       (name2 == "CZ3" && name1 == "CE3"))
                order = 2;
              else if ((name1 == "CZ2" && name2 == "CH2") ||
                       (name2 == "CZ2" && name1 == "CH2"))
                order = 2;
              else if ((name1 == "CE2" && name2 == "CD2") ||
                       (name2 == "CE2" && name1 == "CD2"))
                order = 2;
              break;
            }
            break;
          case 0:
            if ((name1 == "O2P" || name1 == "OP2")) {
              ai1->formalCharge = -1;
              ai1->chemFlag = false;
            } else if ((name2 == "O2P" || name2 == "OP2")) {
              ai2->formalCharge = -1;
              ai2->chemFlag = false;
            }

            if ((name1 == "C2" && name2 == "O2") ||
                (name2 == "C2" && name1 == "O2"))
              order = 2;
            else if ((name1 == "C4" && name2 == "O4") ||
                     (name2 == "C4" && name1 == "O4"))
              order = 2;

            else if ((name1 == "C5" && name2 == "C6") ||
                     (name2 == "C5" && name1 == "C6"))
              order = 2;
            else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                     (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
              order = 2;
            break;
        }
        break;
      case 'U':
        if (resn1[1] == 0) {
          if ((name1 == "O2P" || name1 == "OP2")) {
            ai1->formalCharge = -1;
            ai1->chemFlag = false;
          } else if ((name2 == "O2P" || name2 == "OP2")) {
            ai2->formalCharge = -1;
            ai2->chemFlag = false;
          }
          if ((name1 == "C2" && name2 == "O2") ||
              (name2 == "C2" && name1 == "O2"))
            order = 2;
          else if ((name1 == "C4" && name2 == "O4") ||
                   (name2 == "C4" && name1 == "O4"))
            order = 2;

          else if ((name1 == "C5" && name2 == "C6") ||
                   (name2 == "C5" && name1 == "C6"))
            order = 2;
          else if ((name1 == "P" && (name2 == "O1P" || name2 == "OP1")) ||
                   (name2 == "P" && (name1 == "O1P" || name1 == "OP1")))
            order = 2;
        }
        break;
    }
  }
  *(bond_order) = order;
}

static const char* check_next_pdb_object(const char* p, int skip_to_next) {
  const char* start = p;
  while (*p) {
    if (strstartswith(p, "HEADER")) {
      if (skip_to_next) return p;
      return start;
    } else if (strstartswith(p, "ATOM ") || strstartswith(p, "HETATM")) {
      return start;
    } else if (skip_to_next && strcmp("END", p) == 0) {
      /* if we pass over the END of the current PDB file, then reset start */
      start = p;
    }
    p = nextline(p);
  }
  return nullptr;
}

static int get_multi_object_status(const char* p) { /* expensive -- only call
                                                       this if there is no
                                                       other way to determine
                                                       this information */
  int seen_header = 0;
  while (*p) {
    if (strstartswith(p, "HEADER")) {
      if (seen_header) {
        return 1;
      } else {
        seen_header = true;
      }
    }
    p = nextline(p);
  }
  return -1;
}

/**
 * If any atom in I->AtomInfo contains the wildcard character (from
 * "atom_name_wildcard" or "wildcard" setting), then set the object-level
 * "atom_name_wildcard" setting to " " (disables wildcard matching).
 */
int ObjectMoleculeAutoDisableAtomNameWildcard(ObjectMolecule* I) {
  PyMOLGlobals* G = I->G;
  char wildcard = 0;
  int found_wildcard = false;

  {
    const char* tmp =
        SettingGet_s(G, nullptr, I->Setting.get(), cSetting_atom_name_wildcard);
    if (tmp && tmp[0]) {
      wildcard = *tmp;
    } else {
      tmp = SettingGet_s(G, nullptr, I->Setting.get(), cSetting_wildcard);
      if (tmp) {
        wildcard = *tmp;
      }
    }
    if (wildcard == 32) wildcard = 0;
  }

  if (wildcard) {
    int a;
    const char* p;
    char ch;
    const AtomInfoType* ai = I->AtomInfo;

    for (a = 0; a < I->NAtom; a++) {
      p = LexStr(G, ai->name);
      while ((ch = *(p++))) {
        if (ch == wildcard) {
          found_wildcard = true;
          break;
        }
      }
      ai++;
    }
    if (found_wildcard) {
      ExecutiveSetObjSettingFromString(G, cSetting_atom_name_wildcard, " ", I,
                                       -1, true, true);
    }
  }
  return found_wildcard;
}
#define PDB_MAX_TAGS 64

static void ObjectMoleculePDBStr2CoordSetPASS1(
    PyMOLGlobals* G, int* ok, const char** restart_model, const char* p,
    int n_tags, const char* const* tag_start, int* nAtom, char cc[], int quiet,
    int* bogus_name_alignment, int* ssFlag, const char** next_pdb,
    PDBInfoRec* info, int only_read_one_model, int* ignore_conect,
    int* bondFlag, int* have_bond_order) {
  int seen_end_of_atoms = false;
  *restart_model = nullptr;
  while (*ok && *p) {
    AddOrthoOutputIfMatchesTags(G, n_tags, *nAtom, tag_start, p, cc, quiet);
    if ((strstartswith(p, "ATOM ") || strstartswith(p, "HETATM")) &&
        (!*restart_model)) {
      if (!seen_end_of_atoms) (*nAtom)++;
      if (*bogus_name_alignment) {
        ncopy(cc, nskip(p, 12), 4); /* copy the atom field */
        if ((cc[0] == 32) &&
            (cc[1] != 32)) { /* check to see if indentation was followed
                                correctly, 32 - space */
          *bogus_name_alignment = false;
        }
      }
    } else if (strstartswith(p, "HELIX ")) {
      *ssFlag = true;
    } else if (strstartswith(p, "SHEET ")) {
      *ssFlag = true;
    } else if (strstartswith(p, "HEADER")) {
      if (*nAtom > 0) { /* if we've already found atom records, then this must
                           be a new pdb */
        (*next_pdb) = p;
        break;
      }
    } else if (strstartswith(p, "REMARK")) {
      //        char * start = p;  // USED TO SAVE REMARKS (TODO)
      do {
        if (info && strncmp("    GENERATED BY TRJCONV", p + 6, 24) == 0)
          info->ignore_header_names = true;
        p = nextline(p);
        AddOrthoOutputIfMatchesTags(G, n_tags, *nAtom, tag_start, p, cc, quiet);
      } while (strstartswith(p, "REMARK"));
      // REMARKS are string from start to p, but not saved currently (TODO)
      continue;
    } else if (strstartswith(p, "ENDMDL") && (!*restart_model)) {
      *restart_model = nextline(p);
      seen_end_of_atoms = true;
      if (only_read_one_model) break;
    } else if (strstartswith(p, "END")) { /* stop parsing after END */
      ntrim(cc, p, 6);
      if (strcmp("END", cc) == 0) { /* END */
        seen_end_of_atoms = true;
        if (next_pdb) {
          p = nextline(p);
          ncopy(cc, p, 6);
          if (strcmp("HEADER", cc) == 0) {
            (*next_pdb) = p; /* found another PDB file after this one... */
          } else if (strcmp("ENDMDL", cc) == 0) {
            seen_end_of_atoms = false;
          }
        }
        break;
      }
    } else if (strstartswith(p, "CONECT")) {
      if (!*ignore_conect) *bondFlag = true;
    } else if (strstartswith(p, "USER") && (!*restart_model)) {
    }
    p = nextline(p);
  }
}

/**
 * Datastructure for efficient array-based secondary structure lookup.
 */
struct SSHash {
  int n_ss;          // number of ss_list items
  int* ss[256];      // one array for each chain identifier
  SSEntry* ss_list;  // VLA
};

static void sshash_free(SSHash* hash) {
  int a;
  if (!hash) return;
  for (a = 0; a <= 255; a++) FreeP(hash->ss[a]);
  VLAFreeP(hash->ss_list);
  FreeP(hash);
}

static SSHash* sshash_new() {
  SSHash* hash = pymol::calloc<SSHash>(1);
  ok_assert(1, hash);
  hash->n_ss = 1;
  hash->ss_list = VLAlloc(SSEntry, 50);
  ok_assert(1, hash->ss_list);
  return hash;
ok_except1:
  sshash_free(hash);
  return nullptr;
}

/**
 * Insert a secondary structure record into the hash table.
 */
static int sshash_register_rec(SSHash* hash, unsigned char ss_chain1,
                               int ss_resv1, char ss_inscode1,
                               unsigned char ss_chain2, int ss_resv2,
                               char ss_inscode2, char SSCode) {
  /* pretty confusing how this works... the following efficient (i.e.
     array-based) secondary structure lookup even when there are multiple
     insertion codes and where there may be multiple SS records for the residue
     using different insertion codes */

  unsigned char chain;
  int ss_found = false, ssi = 0, a, b, index;
  SSEntry* sst;

  // up to two iterations:
  // 1) assume chain1==chain2
  // 2) chains are not the same (undefined in PDB spec?)
  for (a = 0, chain = ss_chain1; a < 2; a++, chain = ss_chain2) {
    // allocate new array for chain if necc.
    if (!hash->ss[chain]) {
      ok_assert(1, hash->ss[chain] = pymol::calloc<int>(cResvMask + 1));
    }

    sst = nullptr;
    // iterate over all residues indicated
    for (b = ss_resv1; b <= ss_resv2; b++) {
      index = b & cResvMask;

      if (hash->ss[chain][index]) {
        // make a unique copy in the event of multiple entries for one resv
        sst = nullptr;
      }

      if (!sst) {
        VLACheck(hash->ss_list, SSEntry, hash->n_ss);
        ok_assert(1, hash->ss_list);
        ssi = (hash->n_ss)++;
        sst = hash->ss_list + ssi;
        sst->resv1 = ss_resv1;
        sst->resv2 = ss_resv2;
        sst->chain1 = ss_chain1;
        sst->chain2 = ss_chain2;
        sst->type = SSCode;
        sst->inscode1 = ss_inscode1;
        sst->inscode2 = ss_inscode2;
        ss_found = true;
      }
      sst->next = hash->ss[chain][index];
      hash->ss[chain][index] = ssi;
      if (sst->next) sst = nullptr; /* force another unique copy */
    }
  }
  return ss_found;
ok_except1:
  return false;
}

/**
 * Assign ai->ssType
 */
static void sshash_lookup(SSHash* hash, AtomInfoType* ai,
                          unsigned char ss_chain1) {
  int index, ssi;
  SSEntry* sst = nullptr;

  index = ai->resv & cResvMask;
  if (hash->ss[ss_chain1]) {
    ssi = hash->ss[ss_chain1][index];
    while (ssi) {
      sst = hash->ss_list + ssi;
      /* contains shared entry, or unique linked list for each residue */
      if (ai->resv >= sst->resv1 && ai->resv <= sst->resv2 &&
          (ai->resv != sst->resv1 || ai->inscode >= sst->inscode1) &&
          (ai->resv != sst->resv2 || ai->inscode <= sst->inscode2)) {
        ai->ssType[0] = sst->type;
        return;
      }
      ssi = sst->next;
    }
  }
}

/**
 * PQR atom line parsing
 *
 * Try to parse columns white space delimited (10 columns with optional
 * chain, 9 without).
 *
 * Where PQR files come from:
 *
 * pdb2pqr -> writes PDB-like fixed colums. APBS will fail to read those files
 * if columns are too wide.
 *
 * pdb2pqr --whitespace -> puts extra whitespace, starting at column 2, but
 * not between chain and resi! PyMOL <= 2.1 fails to read those files.
 *
 * APBS Tools Plugin by Michael Lerner adds extra whitespace before negative
 * coordinates (assuming -100.0 is the most likely source of error).
 *
 * @param[in,out] p     Pointer to parse from, points after the ATOM field.
 *                      Will move the pointer to the end of the line if
 *                      parsing was successful.
 * @param[out]    ai    Atom to populate with data
 * @param[out]    coord Coordinates to populate
 * @return true on success, false otherwise.
 */
static bool parse_pqr_atom_line(PyMOLGlobals* G, const char*& p,
                                AtomInfoType* ai, float* coord) {
  auto p_eol = nskip(p, 999);  // end of line pointer
  std::string cc(p, p_eol);    // line (starting after ATOM field)

  // whitespace splitting
  auto columns = strsplit(cc);

  // insert chain column if missing
  if (columns.size() == 9) {
    columns.insert(columns.begin() + 3, "");

    // check for concatenated chain + resi
    if (columns[4].size() > 4 && !isdigit(columns[4][0])) {
      columns[3] = columns[4].substr(0, 1);
      columns[4] = columns[4].substr(1);
    }
  }

  // for validation: if number parsing consumes the entire string, then dummy
  // will never be populated (sscanf(...) == 1, not 2)
  char dummy[2];

  // validate numeric fields and populate atom info and coordinates
  if (columns.size() == 10 &&
      sscanf(columns[0].c_str(), "%d%1s", &ai->id, dummy) == 1 &&
      sscanf(columns[5].c_str(), "%f%1s", coord + 0, dummy) == 1 &&
      sscanf(columns[6].c_str(), "%f%1s", coord + 1, dummy) == 1 &&
      sscanf(columns[7].c_str(), "%f%1s", coord + 2, dummy) == 1 &&
      sscanf(columns[8].c_str(), "%f%1s", &ai->partialCharge, dummy) == 1 &&
      sscanf(columns[9].c_str(), "%f%1s", &ai->elec_radius, dummy) == 1) {
    LexAssign(G, ai->name, columns[1].c_str());
    LexAssign(G, ai->resn, columns[2].c_str());
    LexAssign(G, ai->chain, columns[3].c_str());
    ai->setResi(columns[4].c_str());

    // move parser to next line
    p = p_eol;

    return true;
  }

  return false;
}

CoordSet* ObjectMoleculePDBStr2CoordSet(PyMOLGlobals* G, const char* buffer,
                                        AtomInfoType** atInfoPtr,
                                        const char** restart_model,
                                        char* segi_override, char* pdb_name,
                                        const char** next_pdb, PDBInfoRec* info,
                                        int quiet, int* model_number) {
  const char* p;
  int nAtom;
  int a;
  float* coord = nullptr;
  CoordSet* cset = nullptr;
  AtomInfoType *atInfo = nullptr, *ai;
  int AFlag;
  char SSCode;
  int atomCount;
  int bondFlag = false;
  BondType *bond = nullptr, *ii1, *ii2;
  int* idx;
  int nBond = 0;
  int b1, b2, nReal, maxAt;
  std::unique_ptr<CSymmetry> symmetry;
  int auto_show = RepGetAutoShowMask(G);
  int reformat_names = SettingGetGlobal_i(G, cSetting_pdb_reformat_names_mode);
  int truncate_resn = SettingGetGlobal_b(G, cSetting_pdb_truncate_residue_name);
  const char *tags_in = SettingGetGlobal_s(G, cSetting_pdb_echo_tags),
             *tag_start[PDB_MAX_TAGS];
  int n_tags = 0;
  int foundNextModelFlag = false;
  int ssFlag = false;
  int ss_resv1 = 0, ss_resv2 = 0;
  char ss_inscode1 = '\0', ss_inscode2 = '\0';
  unsigned char ss_chain1 = 0, ss_chain2 = 0;
  SSHash* ss_hash = nullptr;
  char cc[MAXLINELEN], tags[MAXLINELEN];
  int ignore_pdb_segi = 0;
  int ss_valid, ss_found = false;
  int only_read_one_model = false;
  int ignore_conect = SettingGetGlobal_b(G, cSetting_pdb_ignore_conect);
  int have_bond_order = false;
  int seen_model, in_model = false;
  int is_end_of_object = false;
  int literal_names = SettingGetGlobal_b(G, cSetting_pdb_literal_names);
  int bogus_name_alignment = true;
  AtomName literal_name = "";
  int ok = true;
  lexidx_t segi_override_idx = LexIdx(G, segi_override);

  if (tags_in && (!quiet) && (!*restart_model)) {
    char* p = tags;
    strcpy(tags, tags_in);

    while (*p) {
      while (*p == ' ') /* skip spaces */
        p++;
      if (*p) {
        tag_start[n_tags] = p;
        n_tags++;
        while (*p) {
          if (*p != ',') {
            if (*p == ' ') *p = 0;
            p++;
          } else
            break;
        }
        if (*p) { /* terminate tag */
          *p = 0;
          p++;
        }
      }
    }
  }

  if (literal_names) reformat_names = 0;

  ignore_pdb_segi = SettingGetGlobal_b(G, cSetting_ignore_pdb_segi);

  p = buffer;
  nAtom = 0;
  if (atInfoPtr) atInfo = *atInfoPtr;

  if (!atInfo)
    ErrFatal(G, "PDBStr2CoordSet",
             "need atom information record!"); /* failsafe for old version.. */

  if (buffer == *restart_model)
    only_read_one_model = true;
  else if (info && (info->multiplex > 0)) {
    only_read_one_model = true;
    if (!info->multi_object_status) { /* figure out if this is a multi-object
                                         (multi-HEADER) pdb file */
      info->multi_object_status = get_multi_object_status(p);
    }
  }
  /* PASS 1 */
  ObjectMoleculePDBStr2CoordSetPASS1(
      G, &ok, restart_model, p, n_tags, tag_start, &nAtom, cc, quiet,
      &bogus_name_alignment, &ssFlag, next_pdb, info, only_read_one_model,
      &ignore_conect, &bondFlag, &have_bond_order);
  /* INPUT USED:
     only_read_one_model - only reads one pdb model if set
     n_tags, tag_start - tags defined to report in log
     cc - just temporary char * buffer that is used, no input/output
     OUTPUT USED:
     bogus_name_alignment - whether all ATOM/HETATM has correct names in PDB
     (12-16, starting with space ssFlag - if PDB has HELIX and/or SHEET records
     both INPUT/OUTPUT:
     nAtom - returns the number of atoms read in, also uses it to detect
     multiple PDBs ignore_conect - do not set bondFlag if set

     END PASS 1 */

  *restart_model = nullptr;
  if (ok) {
    coord = VLAlloc(float, 3 * nAtom);
    CHECKOK(ok, coord);
  }

  if (ok && atInfo) {
    VLACheck(atInfo, AtomInfoType, nAtom);
    CHECKOK(ok, atInfo);
    if (ok) *atInfoPtr = atInfo;
  }
  if (ok && bondFlag) {
    nBond = 0;
    bond = VLACalloc(BondType, 6 * nAtom);
    CHECKOK(ok, bond);
  }
  p = buffer;
  PRINTFB(G, FB_ObjectMolecule, FB_Blather)
  " ObjectMoleculeReadPDB: Found %i atoms...\n", nAtom ENDFB(G);

  if (ok && ssFlag) {
    ss_hash = sshash_new();
  }

  a = 0; /* WATCHOUT */
  atomCount = 0;

  /* PASS 2 */
  seen_model = false;

  while (ok && *p) {
    /*      printf("%c%c%c%c%c%c\n",p[0],p[1],p[2],p[3],p[4],p[5]); */
    AFlag = false;
    SSCode = 0;
    if (strstartswith(p, "ATOM "))
      AFlag = 1;
    else if (strstartswith(p, "HETATM"))
      AFlag = 2;
    else if (strstartswith(p, "HEADER")) {
      if (pdb_name) {
        if (atomCount > 0) {
          /* have we already read atoms?  then this is probably a new PDB file
           */
          (*next_pdb) = p; /* found another PDB file after this one... */
          break;
        } else if ((!info) || (!info->ignore_header_names)) {
          /* if this isn't an MD trajectory... */
          const char* pp;
          pp = nskip(p, 62); /* is there a four-letter PDB code? */
          pp = ntrim(pdb_name, pp, 4);
          if (pdb_name[0] ==
              0) { /* if not, is there a plain name (for MERCK!) */
            p = nskip(p, 6);
            p = ntrim(cc, p, 44);
            UtilNCopy(pdb_name, cc, WordLength);
          } else {
            p = pp;
          }
        }
      }
    } else if (strstartswith(p, "MODEL ")) {
      if (model_number) {
        int tmp;
        p = nskip(p, 10);
        p = ncopy(cc, p, 5);
        if (sscanf(cc, "%d", &tmp) == 1) *model_number = tmp;
      }
      seen_model = true;
      in_model = true;
    } else if (strstartswith(p, "HELIX ")) {
      ss_valid = true;

      p = nskip(p, 19);
      ss_chain1 = (*p);
      p = nskip(p, 2);
      p = ncopy(cc, p, 4);
      if (!sscanf(cc, "%d", &ss_resv1)) ss_valid = false;
      ss_inscode1 = makeInscode(*p);

      p = nskip(p, 6);
      ss_chain2 = (*p);
      p = nskip(p, 2);
      p = ncopy(cc, p, 4);

      if (!sscanf(cc, "%d", &ss_resv2)) ss_valid = false;
      ss_inscode2 = makeInscode(*p);

      if (ss_valid) {
        PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " ObjectMolecule: read HELIX %c %d%.1s %c %d%.1s\n", ss_chain1,
            ss_resv1, &ss_inscode1, ss_chain2, ss_resv2, &ss_inscode2 ENDFB(G);
        SSCode = 'H';
      }

      if (ss_chain1 == ' ') ss_chain1 = 0;
      if (ss_chain2 == ' ') ss_chain2 = 0;

    } else if (strstartswith(p, "SHEET ")) {
      ss_valid = true;
      p = nskip(p, 21);
      ss_chain1 = (*p);
      p = nskip(p, 1);
      p = ncopy(cc, p, 4);
      if (!sscanf(cc, "%d", &ss_resv1)) ss_valid = false;
      ss_inscode1 = makeInscode(*p);
      p = nskip(p, 6);
      ss_chain2 = (*p);
      p = nskip(p, 1);
      p = ncopy(cc, p, 4);
      if (!sscanf(cc, "%d", &ss_resv2)) ss_valid = false;
      ss_inscode2 = makeInscode(*p);

      if (ss_valid) {
        PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " ObjectMolecule: read SHEET %c %d%.1s %c %d%.1s\n", ss_chain1,
            ss_resv1, &ss_inscode1, ss_chain2, ss_resv2, &ss_inscode2 ENDFB(G);
        SSCode = 'S';
      }

      if (ss_chain1 == ' ') ss_chain1 = 0;
      if (ss_chain2 == ' ') ss_chain2 = 0;

    } else if (strstartswith(p, "ENDMDL")) {
      if (*restart_model)
        in_model = false;
      else {
        *restart_model = nextline(p);
        in_model = false;
        if (only_read_one_model) {
          const char* pp;
          pp = nextline(p);
          if (strstartswith(
                  pp,
                  "END")) { /* END we're going to be starting a new object... */
            (*next_pdb) = check_next_pdb_object(nextline(pp), true);
            if (info && (info->multiplex == 0)) {
              /* multiplex == 0:  FORCED multimodel behavior with concatenated
               * PDB files */
              (*restart_model) = (*next_pdb);
              (*next_pdb) = nullptr;
              foundNextModelFlag = true;
              info->multi_object_status = -1;
            } else {
              is_end_of_object = true;
            }
          } else if (strstartswith(pp, "MODEL")) { /* not a new object...just a
                                                      new state (model) */
            if (info &&
                (info->multiplex > 0)) { /* end object if we're multiplexing */
              (*next_pdb) = check_next_pdb_object(pp, true);
              (*restart_model) = nullptr;
            } else
              is_end_of_object = false;
          } else {
            if (pp[0] > 32) /* more content follows... */
              (*next_pdb) = check_next_pdb_object(pp, true);
            else
              (*next_pdb) = nullptr; /* at end of file */
          }
          break;
        }
      }
    } else if (strstartswith(p, "END")) {
      ntrim(cc, p, 6);
      if ((strcmp("END", cc) == 0) && (!in_model)) { /* stop parsing here... */
        const char* pp;
        pp = nextline(
            p); /* ...unless we're in MODEL or next line is itself ENDMDL */
        p = ncopy(cc, p, 6);
        if (!((cc[0] == 'E') && (cc[1] == 'N') && (cc[2] == 'D') &&
              (cc[3] == 'M') && (cc[4] == 'D') &&
              (cc[5] == 'L'))) { /* NOTE: this test seems unnecessary given
                                    strcmp above... */
          if (!*next_pdb) {
            (*next_pdb) = check_next_pdb_object(pp, false);
          }
          if ((*next_pdb) && info && (!info->multiplex) && !(*restart_model)) {
            /* multiplex == 0:  FORCED multimodel behavior with concatenated PDB
             * files */
            (*restart_model) = (*next_pdb);
            (*next_pdb) = nullptr;
            foundNextModelFlag = true;
            info->multi_object_status = -1;
            is_end_of_object = false;
            break;
          }
          if (*next_pdb) { /* we've found another object... */
            if (*restart_model)
              is_end_of_object =
                  false; /* however, if we're parsing multi-models, then we're
                            not yet at the end */
            else
              is_end_of_object = true;
            break;
          } else if (!seen_model)
            break;
        }
      }
    } else if (strstartswith(p, "CRYST1") && (!*restart_model)) {
      if (!symmetry) {
        symmetry.reset(new CSymmetry(G));
        CHECKOK(ok, symmetry);
      }
      if (symmetry) {
        int symFlag = true;
        PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " PDBStrToCoordSet: Attempting to read symmetry information\n" ENDFB(G);
        p = nskip(p, 6);
        symFlag = true;
        float cellparams[3];
        p = ncopy(cc, p, 9);
        if (sscanf(cc, "%f", cellparams + 0) != 1) symFlag = false;
        p = ncopy(cc, p, 9);
        if (sscanf(cc, "%f", cellparams + 1) != 1) symFlag = false;
        p = ncopy(cc, p, 9);
        if (sscanf(cc, "%f", cellparams + 2) != 1) symFlag = false;

        symmetry->Crystal.setDims(cellparams);

        p = ncopy(cc, p, 7);
        if (sscanf(cc, "%f", cellparams + 0) != 1) symFlag = false;
        p = ncopy(cc, p, 7);
        if (sscanf(cc, "%f", cellparams + 1) != 1) symFlag = false;
        p = ncopy(cc, p, 7);
        if (sscanf(cc, "%f", cellparams + 2) != 1) symFlag = false;

        symmetry->Crystal.setAngles(cellparams);

        p = nskip(p, 1);
        p = ncopy(cc, p, 11);
        UtilCleanStr(cc);
        symmetry->setSpaceGroup(cc);
        p = ncopy(cc, p, 3);
        if (sscanf(cc, "%d", &symmetry->PDBZValue) != 1)
          symmetry->PDBZValue = 1;
        if (!symFlag) {
          ErrMessage(G, "PDBStrToCoordSet", "Error reading CRYST1 record\n");
          symmetry.reset();
        }
      }
    } else if (strstartswith(p, "SCALE") && (!*restart_model) &&
               info) { /* SCALEn */
      switch (p[5]) {
        case '1':
        case '2':
        case '3': {
          int flag = (p[5] - '1');
          int offset = flag * 4;
          int scale_flag = true;
          p = nskip(p, 10);
          p = ncopy(cc, p, 10);
          if (sscanf(cc, "%f", &info->scale.matrix[offset]) != 1)
            scale_flag = false;
          p = ncopy(cc, p, 10);
          if (sscanf(cc, "%f", &info->scale.matrix[offset + 1]) != 1)
            scale_flag = false;
          p = ncopy(cc, p, 10);
          if (sscanf(cc, "%f", &info->scale.matrix[offset + 2]) != 1)
            scale_flag = false;
          p = nskip(p, 5);
          p = ncopy(cc, p, 10);
          if (sscanf(cc, "%f", &info->scale.matrix[offset + 3]) != 1)
            scale_flag = false;
          if (scale_flag) info->scale.flag[flag] = true;
          PRINTFB(G, FB_ObjectMolecule, FB_Blather)
          " PDBStrToCoordSet: SCALE%d %8.4f %8.4f %8.4f %8.4f\n", flag + 1,
              info->scale.matrix[offset], info->scale.matrix[offset + 1],
              info->scale.matrix[offset + 2],
              info->scale.matrix[offset + 3] ENDFB(G);
        } break;
      }
    } else if (strstartswith(p, "CONECT") && bondFlag && (!ignore_conect) &&
               ((!*restart_model) || (!in_model))) {
      p = nskip(p, 6);
      p = ncopy(cc, p, 5);
      if (sscanf(cc, "%d", &b1) == 1)
        while (1) {
          p = ncopy(cc, p, 5);
          if (sscanf(cc, "%d", &b2) != 1)
            break;
          else {
            if ((b1 >= 0) && (b2 >= 0) &&
                (b1 != b2)) { /* IDs must be positive and distinct */
              VLACheck(bond, BondType, nBond);
              CHECKOK(ok, bond);
              if (ok) {
                if (b1 <= b2) {
                  bond[nBond].index[0] =
                      b1; /* temporarily store the atom indexes */
                  bond[nBond].index[1] = b2;
                  bond[nBond].order = 1;
                } else {
                  bond[nBond].index[0] = b2;
                  bond[nBond].index[1] = b1;
                  bond[nBond].order = 1;
                }
                nBond++;
              }
            }
          }
        }
    } else if (strstartswith(p, "USER") && (!*restart_model)) {
    } else if (strstartswith(p, "ANISOU") && (!*restart_model) && (atomCount)) {
      ai = atInfo + atomCount - 1;

      /* TODO: check atom identifier match */

      {
        int dummy;
        p = nskip(p, 6);
        p = ncopy(cc, p, 5);
        if (!sscanf(cc, "%d", &dummy)) dummy = 0;
        if (dummy == ai->id) { /* ATOM ID must match */
          int dummy;
          float* anisou = ai->get_anisou();
          p = nskip(p, 17);
          for (int i = 0; i < 6; ++i) {
            p = ncopy(cc, p, 7);
            if (sscanf(cc, "%d", &dummy)) anisou[i] = dummy / 10000.0F;
          }
        }
      }
    }

    /* END KEYWORDS */

    /* Secondary structure records */

    if (ok && SSCode) {
      ss_found = sshash_register_rec(ss_hash, ss_chain1, ss_resv1, ss_inscode1,
                                     ss_chain2, ss_resv2, ss_inscode2, SSCode);
    }
    /* Atom records */

    if (ok && AFlag && (!*restart_model)) {
      ai = atInfo + atomCount;
      p = nskip(p, 6);

      ai->rank = atomCount;

      if (info && info->is_pqr_file()) {
        if (parse_pqr_atom_line(G, p, ai, coord + a)) {
          goto pqr_done;
        }
      }

      p = ncopy(cc, p, 5);
      if (!sscanf(cc, "%d", &ai->id)) ai->id = 0;

      p = nskip(p, 1); /* to 12 */
      p = ncopy(literal_name, p, 4);
      if (literal_names) {
        LexAssign(G, ai->name, literal_name);
      } else {
        ParseNTrim(cc, literal_name, 4);
        LexAssign(G, ai->name, cc);
      }

      p = ncopy(cc, p, 1);
      if (*cc == 32)
        ai->alt[0] = 0;
      else {
        ai->alt[0] = *cc;
        ai->alt[1] = 0;
      }

      p = ntrim(cc, p, 4); /* now allowing for 4-letter residues */
      if (truncate_resn)   /* unless specifically disabled */
        cc[3] = 0;

      LexAssign(G, ai->resn, cc);

      if (ai->name) {
        const char* ai_name = LexStr(G, ai->name);
        int name_len = strlen(ai_name);
        char name[5];
        switch (reformat_names) {
          case 1: /* pdb compliant: HH12 becomes 2HH1, etc. */
            if (name_len > 3) {
              if ((ai_name[0] >= 'A') && ((ai_name[0] <= 'Z')) &&
                  isdigit(ai_name[3])) {
                if (!(((ai_name[1] >= 'a') && (ai_name[1] <= 'z')) ||
                      ((ai_name[0] == 'C') &&
                       (ai_name[1] == 'L')) || /* try to be smart about */
                      ((ai_name[0] == 'B') &&
                       (ai_name[1] == 'R')) || /* distinguishing common atoms */
                      ((ai_name[0] == 'C') &&
                       (ai_name[1] == 'A')) || /* in all-caps from typical */
                      ((ai_name[0] == 'F') &&
                       (ai_name[1] == 'E')) || /* nonatomic abbreviations */
                      ((ai_name[0] == 'C') && (ai_name[1] == 'U')) ||
                      ((ai_name[0] == 'N') && (ai_name[1] == 'A')) ||
                      ((ai_name[0] == 'N') && (ai_name[1] == 'I')) ||
                      ((ai_name[0] == 'M') && (ai_name[1] == 'G')) ||
                      ((ai_name[0] == 'M') && (ai_name[1] == 'N')) ||
                      ((ai_name[0] == 'H') && (ai_name[1] == 'G')) ||
                      ((ai_name[0] == 'S') && (ai_name[1] == 'E')) ||
                      ((ai_name[0] == 'S') && (ai_name[1] == 'I')) ||
                      ((ai_name[0] == 'Z') && (ai_name[1] == 'N')))) {
                  strncpy(name + 1, ai_name, 3);
                  name[0] = ai_name[3];
                  name[4] = 0;
                  LexAssign(G, ai->name, name);
                }
              }
            } else if (name_len == 3) {
              if ((ai_name[0] == 'H') && (ai_name[1] >= 'A') &&
                  ((ai_name[1] <= 'Z')) && isdigit(ai_name[2])) {
                AtomInfoGetPDB3LetHydroName(G, LexStr(G, ai->resn), ai_name,
                                            name);
                LexAssign(G, ai->name, (name[0] == ' ') ? (name + 1) : name);
              }
            }
            break;
          case 2: /* amber compliant: 2HH1 becomes HH12 */
          case 3: /* pdb compliant, but use IUPAC within PyMOL */
            if (ai_name[0]) {
              if (isdigit(ai_name[0]) && ai_name[1] && (!isdigit(ai_name[1]))) {
                if (1 < name_len && name_len < 5) {
                  strcpy(name, ai_name + 1);
                  name[name_len - 1] = ai_name[0];
                  name[name_len] = 0;
                  LexAssign(G, ai->name, name);
                }
                break;
                default: /* AS IS */
                  break;
              }
            }
            break;
          case 4: /* simply read trim and write back out with 3-letter names
                     starting from the second column, and four-letter names
                     starting in the first */
            ntrim(cc, ai_name, 4);
            LexAssign(G, ai->name, cc);
            break;
        }
      }

      p = ncopy(cc, p, 1);
      if (ai->chain) {
        LexDec(G, ai->chain);
      }
      if (*cc == ' ') {
        ss_chain1 = 0;
        ai->chain = 0;
      } else {
        ss_chain1 = *cc;
        ai->chain = LexIdx(G, cc);
      }

      p = ncopy(cc, p, 4);
      if (!sscanf(cc, "%d", &ai->resv)) ai->resv = 0;
      ai->setInscode(*p);
      p = nskip(p, 1);

      if (ssFlag) { /* get secondary structure information (if avail) */
        sshash_lookup(ss_hash, ai, ss_chain1);
      } else {
        ai->cartoon = cCartoon_tube;
      }

      {
        p = nskip(p, 3);
        p = ncopy(cc, p, 8);
        sscanf(cc, "%f", coord + a);
        p = ncopy(cc, p, 8);
        sscanf(cc, "%f", coord + (a + 1));
        p = ncopy(cc, p, 8);
        sscanf(cc, "%f", coord + (a + 2));
      }

      if ((!info) || (!info->is_pqr_file())) { /* standard PDB file */
        p = ncopy(cc, p, 6);
        if (!sscanf(cc, "%f", &ai->q)) ai->q = 1.0;

        p = ncopy(cc, p, 6);
        if (!sscanf(cc, "%f", &ai->b)) ai->b = 0.0;

        if (info->variant == PDB_VARIANT_PDBQT) {
          ignore_pdb_segi = true;
          p = nskip(p, 4);
          p = ncopy(cc, p, 6);
          if (!sscanf(cc, "%f", &ai->partialCharge)) ai->partialCharge = 0.0;

          // type is 78-79 in pdbqt, 77-78 in pdb
          p = nskip(p, 1);
        } else {
          p = nskip(p, 6);
          p = ncopy(cc, p, 4);
        }

        if (!ignore_pdb_segi) {
          if (!segi_override_idx) {
            if (cc[3] == '1' && atomCount && strncmp(p, "0000", 4) == 0) {
              /* atom ID overflow? (nonstandard use...)... */
              LexAssign(G, segi_override_idx, (ai - 1)->segi);
              LexAssign(G, ai->segi, (ai - 1)->segi);
            } else {
              UtilCleanStr(cc);
              LexAssign(G, ai->segi, cc);
            }
          } else {
            LexAssign(G, ai->segi, segi_override_idx);
          }
        } else {
          LexAssign(G, ai->segi, 0);
        }

        p = ncopy(cc, p, 2);
        if (!sscanf(cc, "%s", ai->elem))
          ai->elem[0] = 0;
        else if (!((((ai->elem[0] >= 'a') &&
                     (ai->elem[0] <=
                      'z')) || /* don't get confused by PDB misuse */
                    ((ai->elem[0] >= 'A') && (ai->elem[0] <= 'Z'))) &&
                   (((ai->elem[1] == 0) ||
                     ((ai->elem[1] >= 'a') && (ai->elem[1] <= 'z')) ||
                     ((ai->elem[1] >= 'A') && (ai->elem[1] <= 'Z'))))))
          ai->elem[0] = 0;
        else if (info->variant == PDB_VARIANT_PDBQT) {
          if (strcmp(ai->elem, "A") == 0) {
            // aromatic carbon
            ai->elem[0] = 'C';
          } else if (isupper(ai->elem[1])) {
            // h-bond donor or acceptor
            ai->elem[1] = 0;
          }
        }

        if (!ai->elem[0]) {
          if (((literal_name[0] == ' ') ||
               ((literal_name[0] >= '0') && (literal_name[0] <= '9'))) &&
              (literal_name[1] >= 'A') &&
              (literal_name[1] <= 'Z')) { /* infer element from name column */
            ai->elem[0] = literal_name[1];
            ai->elem[1] = 0;
          } else if (((literal_name[0] >= 'A') && (literal_name[0] <= 'Z')) &&
                     (((literal_name[1] >= 'A') && (literal_name[1] <= 'Z')) ||
                      ((literal_name[1] >= 'a') &&
                       (literal_name[1] <=
                        'z')))) { /* infer element from name column */
            ai->elem[0] = literal_name[0];
            ai->elem[2] = 0;
            if ((literal_name[1] >= 'A') &&
                (literal_name[1] <= 'Z')) { /* second letter is capitalized */
              if (bogus_name_alignment) {
                /* if other atom names aren't properly aligned */
                ai->elem[1] = 0; /* kill 2nd letter */
              } else if (literal_name[0] == 'H') {
                /* or if this is an ultra-bogus PDB with inconsistent
                   indendentation, and this is likely a hydrogen */
                ai->elem[1] = 0; /* kill 2nd letter */
              } else {
                ai->elem[1] = tolower(literal_name[1]);
              }
            } else
              ai->elem[1] = literal_name[1];
          }
        }

        p = ncopy(cc, p, 2);
        if ((cc[1] == '-') || (cc[1] == '+')) {
          /* only read formal charge when sign is present */
          char ctmp = cc[0];
          cc[0] = cc[1];
          cc[1] = ctmp;
          if (!sscanf(cc, "%hhi", &ai->formalCharge)) ai->formalCharge = 0;
        }

        /* end normal PDB */
      } else if (info && info->is_pqr_file()) {
        p = ParseWordNumberCopy(cc, p, MAXLINELEN - 1);
        if (!sscanf(cc, "%f", &ai->partialCharge)) ai->partialCharge = 0.0F;

        p = ParseWordNumberCopy(cc, p, MAXLINELEN - 1);
        if (sscanf(cc, "%f", &ai->elec_radius) != 1) ai->elec_radius = 0.0F;
      }

    pqr_done:

      ai->visRep = auto_show;

      if (AFlag == 1)
        ai->hetatm = 0;
      else {
        ai->hetatm = 1;
        ai->flags = cAtomFlag_ignore;
      }

      AtomInfoAssignParameters(G, ai);
      AtomInfoAssignColors(G, ai);

      PRINTFD(G, FB_ObjectMolecule)
      "%s %s %d%c %s %8.3f %8.3f %8.3f %6.2f %6.2f %s\n", LexStr(G, ai->name),
          LexStr(G, ai->resn), ai->resv, ai->getInscode(true),
          LexStr(G, ai->chain), *(coord + a), *(coord + a + 1),
          *(coord + a + 2), ai->b, ai->q, LexStr(G, ai->segi) ENDFD;

      if (atomCount < nAtom) { /* safety */
        a += 3;
        atomCount++;
      }
    }
    p = nextline(p);
  }

  /* END PASS 2 */

  if (ok && bondFlag) {
    UtilSortInPlace(G, bond, nBond, sizeof(BondType),
                    (UtilOrderFn*)BondInOrder);
    if (nBond) {
      if (!have_bond_order) { /* handle PDB bond-order kludge */
        ii1 = bond;
        ii2 = bond + 1;
        nReal = 1;
        ii1->order = 1;
        a = nBond - 1;
        while (a) {
          if ((ii1->index[0] == ii2->index[0]) &&
              (ii1->index[1] == ii2->index[1])) {
            ii1->order++; /* count dup */
          } else {
            ii1++; /* non-dup, make copy */
            ii1->index[0] = ii2->index[0];
            ii1->index[1] = ii2->index[1];
            ii1->order = ii2->order;
            nReal++;
          }
          ii2++;
          a--;
        }
        nBond = nReal;
      }
      /* now, find atoms we're looking for */

      /* determine ranges */
      maxAt = nAtom;
      ii1 = bond;
      for (a = 0; a < nBond; a++) {
        if (ii1->index[0] > maxAt) maxAt = ii1->index[0];
        if (ii1->index[1] > maxAt) maxAt = ii1->index[1];
        ii1++;
      }
      for (a = 0; a < nAtom; a++)
        if (maxAt < atInfo[a].id) maxAt = atInfo[a].id;
      /* build index */
      maxAt++;
      idx = pymol::malloc<int>(maxAt + 1);
      CHECKOK(ok, idx);
      if (ok) {
        for (a = 0; a < maxAt; a++) {
          idx[a] = -1;
        }
        for (a = 0; a < nAtom; a++) idx[atInfo[a].id] = a;

        /* convert indices to bonds */
        ii1 = bond;
        ii2 = bond;
        nReal = 0;
      }
      if (ok) {
        int unbond_cations = SettingGetGlobal_i(G, cSetting_pdb_unbond_cations);
        int flag;

        for (a = 0; a < nBond; a++) {
          if ((ii1->index[0] >= 0) && ((ii1->index[1]) >= 0)) {
            if ((idx[ii1->index[0]] >= 0) &&
                (idx[ii1->index[1]] >=
                 0)) { /* in case PDB file has bad bonds */
              ii2->index[0] = idx[ii1->index[0]];
              ii2->index[1] = idx[ii1->index[1]];
              ii2->order = ii1->order;
              if ((ii2->index[0] >= 0) && (ii2->index[1] >= 0)) {
                if (!have_bond_order) { /* handle PDB bond order kludge */
                  if (ii1->order <= 2)
                    ii2->order = 1;
                  else if (ii1->order <= 4)
                    ii2->order = 2;
                  else if (ii1->order <= 6)
                    ii2->order = 3;
                  else
                    ii2->order = 4;
                }
                flag = true;
                if (unbond_cations) {
                  if (AtomInfoIsFreeCation(G, atInfo + ii2->index[0]))
                    flag = false;
                  else if (AtomInfoIsFreeCation(G, atInfo + ii2->index[1]))
                    flag = false;
                }
                if (flag) {
                  atInfo[ii2->index[0]].bonded = true;
                  atInfo[ii2->index[1]].bonded = true;
                  nReal++;
                  ii2++;
                }
              }
            }
          }
          ii1++;
        }
      }
      nBond = nReal;
      FreeP(idx);
    }
  }
  if (ss_found && !quiet) {
    PRINTFB(G, FB_ObjectMolecule, FB_Details)
    " ObjectMolecule: Read secondary structure assignments.\n" ENDFB(G);
  }
  if (symmetry && !quiet && (!only_read_one_model)) {
    PRINTFB(G, FB_ObjectMolecule, FB_Details)
    " ObjectMolecule: Read crystal symmetry information.\n" ENDFB(G);
  }
  PRINTFB(G, FB_ObjectMolecule, FB_Blather)
  " PDBStr2CoordSet: Read %d bonds from CONECT records (%p).\n", nBond,
      (void*)bond ENDFB(G);

  if (ok) {
    cset = CoordSetNew(G);
    CHECKOK(ok, cset);

    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->TmpBond = pymol::vla_take_ownership(bond);
    cset->NTmpBond = nBond;
    cset->Symmetry = std::move(symmetry);
    if (atInfoPtr) *atInfoPtr = atInfo;

    if ((*restart_model) && (!foundNextModelFlag)) {
      /* if plan on need to reading another model into this object,
         make sure there is another model to read... */
      p = *restart_model;
      while (*p) {
        if (strstartswith(p, "HEADER")) {
          /* seeing HEADER means we're off the end of the existing file */
          break;
        } else if (strstartswith(p, "MODEL ") || strstartswith(p, "ENDMDL")) {
          foundNextModelFlag = true;
          break;
        }
        p = nextline(p);
      }
      if (!foundNextModelFlag) {
        *restart_model = nullptr;
      }
    }
  }
  if (!ok) {
    if (cset) {
      delete cset;
      cset = nullptr;
    } else {
      VLAFreeP(coord);
      VLAFreeP(bond);
    }
  }
  sshash_free(ss_hash);

  if (ok) {
    if (!seen_model) *model_number = 1;

    if ((*restart_model) &&
        (*next_pdb)) { /* if we're mixing multistate objects and
                          trajectories, then enforce sanity by
                          reading the models first... */
      if (is_end_of_object)
        (*restart_model) = nullptr;
      else if ((*next_pdb) < (*restart_model))
        (*next_pdb) = nullptr;
    }
  }

  LexDec(G, segi_override_idx);
  return (cset);
}
