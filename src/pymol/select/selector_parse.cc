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

// The selection language: keywords and parsing expressions.

#include "pymol/base/feedback.h"
#include "pymol/select/selector.h"
#include "pymol/select/selector_internal.h"

#define STYP_VALU 0
#define STYP_PVAL 0

bool SelectorNameIsKeyword(PyMOLGlobals* G, const char* name) {
  auto I = G->SelectorMgr;
  std::string lowername = name;
  std::transform(lowername.begin(), lowername.end(), lowername.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return I->Key.count(lowername) != 0;
}

/**
 * Removes matching quotes from a string, at string start as well as after word
 * list separators ("+" and ","). Does not consider backslash escaping.
 *
 * Examples (not sure if all of these are intentional):
 * @verbatim
   "foo bar" -> foo bar
   'foo bar' -> foo bar
   "foo"+'bar' -> foo+bar
   "foo bar\" -> foo bar\       # backslash has no escape function
   "foo" "bar" -> foo "bar"     # second pair of quotes not after separator
   foo''+''bar -> foo''+bar     # first pair of quotes not after separator
   "foo"bar" -> foobar"         # third quote unmatched
   foo'+'bar -> foo'+'bar       # no matching quotes after separator
   @endverbatim
 */
static void remove_quotes(std::string& str) {
  /* nasty */

  char* st = &str[0];
  char *p, *q;
  char* quote_start = nullptr;
  char active_quote = 0;
  p = st;
  q = st;

  while (*p) {
    if (((*p) == 34) || ((*p) == 39)) {
      if (quote_start && (active_quote == *p)) { /* eliminate quotes... */
        while (quote_start < (q - 1)) {
          *(quote_start) = *(quote_start + 1);
          quote_start++;
        }
        q--;
        quote_start = nullptr;
        p++;
        continue;
      } else if (quote_start) {
      } else {
        if (p == st) { /* at start => real quote */
          quote_start = q;
          active_quote = *p;
        } else if ((*(p - 1) == '+') ||
                   (*(p - 1) == ',')) { /* after separator => real quote */
          quote_start = q;
          active_quote = *p;
        }
      }
    }
    if (q < p) {
      *q = *p;
    }
    ++q;
    ++p;
  }
  if (q < p) {
    str.resize(q - st);
  }
}

#define STACK_PUSH_VALUE(value)                \
  {                                            \
    depth++;                                   \
    VecCheck(Stack, depth);                    \
    e = Stack.data() + depth;                  \
    e->level = (level << 4) + 1;               \
    e->imp_op_level = (imp_op_level << 4) + 1; \
    imp_op_level = level;                      \
    e->type = STYP_VALU;                       \
    e->m_text = value;                         \
    remove_quotes(e->m_text);                  \
  }

#define STACK_PUSH_OPERATION(ocode)                    \
  {                                                    \
    depth++;                                           \
    VecCheck(Stack, depth);                            \
    e = Stack.data() + depth;                          \
    e->code = ocode;                                   \
    e->level = (level << 4) + ((e->code & 0xF0) >> 4); \
    e->imp_op_level = (imp_op_level << 4) + 1;         \
    imp_op_level = level;                              \
    e->type = (e->code & 0xF);                         \
  }

/**
 * Indicates at wich token the parsing stopped
 * @param tokens token list
 * @param pos token index
 * @return " ".join(tokens[:pos + 1]) + "<--"
 */
static std::string indicate_last_token(const std::vector<std::string>& tokens,
                                       int pos) {
  std::string msg;
  for (int i = 0, i_end = std::min<int>(pos + 1, tokens.size());  //
       i < i_end; ++i) {
    if (i && tokens[i][0]) {
      msg += " ";
    }
    msg += tokens[i];
  }
  msg += "<--";
  return msg;
}

#define return_error_with_tokens(msg) \
  return pymol::make_error(msg, "\n", indicate_last_token(word, c))

#define return_on_error_with_tokens(expr)                    \
  {                                                          \
    auto _temp_result = (expr);                              \
    if (!_temp_result) {                                     \
      return_error_with_tokens(_temp_result.error().what()); \
    }                                                        \
  }
pymol::Result<sele_array_t> SelectorEvaluate(PyMOLGlobals* G,
                                             std::vector<std::string>& word,
                                             int state, int quiet) {
  int level = 0, imp_op_level = 0;
  int depth = 0;
  int a, b, c = 0;
  int ok = true;
  unsigned int code = 0;
  int valueFlag = 0; /* are we expecting? */
  int opFlag, maxLevel;
  int totDepth = 0;
  int exact = 0;

  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  /* CFGs can efficiently be parsed by stacks; use a clean stack w/space
   * for 10 (was: 100) elements */
  EvalElem* e;
  auto Stack = std::vector<EvalElem>(10);

  /* converts all keywords into code, adds them into a operation list */
  while (ok && c < word.size()) {
    if (word[c][0] == '#') {
      if ((!valueFlag) && (!level)) {
        word.resize(c); /* terminate selection if we encounter a comment */
        break;
      }
    }
    switch (word[c][0]) {
      case 0:
        break;
      case '(':
        /* increase stack depth on open parens: (selection ((and token) blah))
         */
        if (valueFlag) return_error_with_tokens("Misplaced (.");
        if (ok) level++;
        break;
      case ')':
        /* decrease stack depth */
        if (valueFlag) return_error_with_tokens("Misplaced ).");
        if (ok) {
          level--;
          if (level < 0)
            return_error_with_tokens("Syntax error.");
          else
            imp_op_level = level;
        }
        if (ok && depth) Stack[depth].level--;
        break;
      default:
        if (valueFlag > 0) { /* standard operand */
          STACK_PUSH_VALUE(word[c]);
          valueFlag--;
        } else if (valueFlag < 0) { /* operation parameter i.e. around X<-- */
          depth++;
          VecCheck(Stack, depth);
          e = Stack.data() + depth;
          e->level = (level << 4) + 1;
          e->imp_op_level = (imp_op_level << 4) + 1;
          imp_op_level = level;
          e->type = STYP_PVAL;
          e->m_text = word[c];
          valueFlag++;
        } else { /* possible keyword... */
          code = WordKey(G, Keyword, word[c].c_str(), 4, ignore_case, &exact);
          if (!code) {
            b = word[c].size() - 1;
            if ((b > 2) && (word[c][b] == ';')) {
              /* kludge to accomodate unnec. ';' usage */
              word[c].resize(b);
              code =
                  WordKey(G, Keyword, word[c].c_str(), 4, ignore_case, &exact);
            } else if (!word[c].compare(0, 2, "p.")) {
              // kludge to parse p.propertyname without space after p.
              code = SELE_PROP;
              exact = 1;
              word[c].erase(0, 2);
              c--;
            }
          }
          PRINTFD(G, FB_Selector)
          " Selector: code %x\n", code ENDFD;
          if ((code > 0) && (!exact))
            if (SelectorIndexByName(G, word[c].c_str()) >= 0)
              code = 0; /* favor selections over partial keyword matches */
          if (code) {
            /* this is a known operation */
            STACK_PUSH_OPERATION(code);
            switch (e->type) {
              case STYP_SEL0:
                valueFlag = 0;
                break;
              case STYP_SEL1:
                valueFlag = 1;
                break;
              case STYP_SEL2:
                valueFlag = 2;
                break;
              case STYP_SEL3:
                valueFlag = 3;
                break;
              case STYP_OPR1:
                valueFlag = 0;
                break;
              case STYP_OPR2:
                valueFlag = 0;
                break;
              case STYP_PRP1:
                valueFlag = -1;
                break;
              case STYP_OP22:
                valueFlag = -2;
                break;
            }
          } else {
            if ((a = std::count(word[c].begin(), word[c].end(),
                                '/'))) { /* handle slash notation */
              if (a > 5) {
                return_error_with_tokens("too many slashes in selection macro");
              }

              // macro codes (some special cases apply! see code below)
              const int macrocodes[] = {SELE_SELs, SELE_SEGs, SELE_CHNs,
                                        SELE_RSNs, SELE_NAMs};

              // two code/value pairs to support resn`resi and name`alt
              int codes[] = {0, 0, 0};  // null-terminated
              char* values[2];

              std::string tmpKW = word[c];
              char* q = &tmpKW[0];

              // if macro starts with "/" then read from left, otherwise
              // read from right
              if (*q == '/') {
                a = 4;
                q++;
              }

              // loop over macro elements
              for (b = 0; q && *q; a--) {
                values[0] = q;

                // null-terminate current element
                if ((q = strchr(q, '/'))) *(q++) = '\0';

                // skip empty elements
                if (!*values[0]) continue;

                codes[0] = macrocodes[4 - a];
                codes[1] = 0;

                // resn`resi or name`alt
                if (codes[0] == SELE_RSNs || codes[0] == SELE_NAMs) {
                  char* backtick = strchr(values[0], '`');
                  if (backtick) {
                    if (codes[0] == SELE_RSNs) {
                      codes[1] = SELE_RSIs;
                    } else {
                      codes[1] = SELE_ALTs;
                    }
                    values[1] = backtick + 1;
                    *backtick = '\0';
                  } else if (codes[0] == SELE_RSNs && values[0][0] >= '0' &&
                             values[0][0] <= '9') {
                    // numeric -> resi
                    codes[0] = SELE_RSIs;
                  }
                }

                for (int i = 0; codes[i]; ++i) {
                  // skip empty and "*" values
                  if (*values[i] && strcmp(values[i], "*")) {
                    if (b++) STACK_PUSH_OPERATION(SELE_AND2);

                    STACK_PUSH_OPERATION(codes[i]);
                    STACK_PUSH_VALUE(values[i]);
                  }
                }
              }

              // all-empty slash macro equals "all"
              if (!b) STACK_PUSH_OPERATION(SELE_ALLz);

            } else if (word[c].find('`') !=
                       std::string::npos) { /* handle <object`index> syntax */
              STACK_PUSH_OPERATION(SELE_MODs);
              valueFlag = 1;
              c--;
            } else { /* handle <selection-name> syntax */
              STACK_PUSH_OPERATION(SELE_SELs);
              valueFlag = 1;
              c--;
            }
          }
        }
        break;
    }
    if (ok) c++; /* go onto next word */
  }
  if (level > 0) {
    return_error_with_tokens("Malformed selection.");
  }
  if (ok) { /* this is the main operation loop */
    totDepth = depth;
    opFlag = true;
    maxLevel = -1;
    for (a = 1; a <= totDepth; a++) {
      if (Stack[a].level > maxLevel) maxLevel = Stack[a].level;
    }
    level = maxLevel;
    PRINTFD(G, FB_Selector)
    " Selector: maxLevel %d %d\n", maxLevel, totDepth ENDFD;
    if (level >= 0)
      while (ok) { /* loop until all ops at all levels have been tried */

        /* order & efficiency of this algorithm could be improved... */

        PRINTFD(G, FB_Selector)
        " Selector: new cycle...\n" ENDFD;
        depth = 1;
        opFlag = true;
        while (ok && opFlag) { /* loop through all entries looking for ops at
                                  the current level */
          opFlag = false;

          if (Stack[depth].level >= level) {
            Stack[depth].level = level; /* trim peaks */
          }
          if (ok)
            if (depth > 0)
              if ((!opFlag) && (Stack[depth].type == STYP_SEL0)) {
                opFlag = true;
                ok = SelectorSelect0(G, &Stack[depth]);
              }
          if (ok)
            if (depth > 1)
              if (Stack[depth - 1].level >= Stack[depth].level) {
                if (ok && (!opFlag) && (Stack[depth - 1].type == STYP_SEL1) &&
                    (Stack[depth].type == STYP_VALU)) {
                  /* 1 argument selection operator */
                  opFlag = true;
                  return_on_error_with_tokens(
                      SelectorSelect1(G, &Stack[depth - 1], quiet));
                  for (a = depth + 1; a <= totDepth; a++)
                    Stack[a - 1] = std::move(Stack[a]);
                  totDepth--;
                } else if (ok && (!opFlag) &&
                           (Stack[depth - 1].type == STYP_OPR1) &&
                           (Stack[depth].type == STYP_LIST)) {
                  /* 1 argument logical operator */
                  opFlag = true;
                  ok = SelectorLogic1(G, &Stack[depth - 1], state);
                  for (a = depth + 1; a <= totDepth; a++)
                    Stack[a - 1] = std::move(Stack[a]);
                  totDepth--;
                } else if ((Stack[depth - 1].type == STYP_LIST) &&
                           (Stack[depth].type == STYP_LIST) &&
                           (!((Stack[depth - 1].level & 0xF) ||
                              (Stack[depth].level & 0xF)))) {
                  /* two adjacent lists at zeroth priority level
                     for the scope (lowest nibble of level is
                     zero) is an implicit OR action */
                  VecCheck(Stack, totDepth + 1);
                  for (a = totDepth; a >= depth; a--)
                    Stack[a + 1] = std::move(Stack[a]);
                  totDepth++;
                  Stack[depth].type = STYP_OPR2;
                  Stack[depth].code = SELE_IOR2;
                  Stack[depth].level = Stack[depth].imp_op_level;
                  Stack[depth].m_text.clear();
                  if (level < Stack[depth].level) level = Stack[depth].level;
                  opFlag = true;
                }
              }
          if (ok)
            if (depth > 2)
              if ((Stack[depth - 1].level >= Stack[depth].level) &&
                  (Stack[depth - 1].level >= Stack[depth - 2].level)) {
                if (ok && (!opFlag) && (Stack[depth - 1].type == STYP_OPR2) &&
                    (Stack[depth].type == STYP_LIST) &&
                    (Stack[depth - 2].type == STYP_LIST)) {
                  /* 2 argument logical operator */
                  ok = SelectorLogic2(G, &Stack[depth - 2]);
                  opFlag = true;
                  for (a = depth + 1; a <= totDepth; a++)
                    Stack[a - 2] = std::move(Stack[a]);
                  totDepth -= 2;
                } else if (ok && (!opFlag) &&
                           (Stack[depth - 1].type == STYP_PRP1) &&
                           (Stack[depth].type == STYP_PVAL) &&
                           (Stack[depth - 2].type == STYP_LIST)) {
                  /* 2 argument logical operator */
                  ok = SelectorModulate1(G, &Stack[depth - 2], state);
                  opFlag = true;
                  for (a = depth + 1; a <= totDepth; a++)
                    Stack[a - 2] = std::move(Stack[a]);
                  totDepth -= 2;
                }
              }
          if (ok)
            if (depth > 2)
              if ((Stack[depth - 2].level >= Stack[depth - 1].level) &&
                  (Stack[depth - 2].level >= Stack[depth].level)) {
                if (ok && (!opFlag) && (Stack[depth - 2].type == STYP_SEL2) &&
                    (Stack[depth - 1].type == STYP_VALU) &&
                    (Stack[depth].type == STYP_VALU)) {
                  /* 2 argument value operator */
                  ok = SelectorSelect2(G, &Stack[depth - 2], state);
                  opFlag = true;
                  for (a = depth + 1; a <= totDepth; a++)
                    Stack[a - 2] = std::move(Stack[a]);
                  totDepth -= 2;
                }
              }
          if (ok)
            if (depth > 3)
              if ((Stack[depth - 3].level >= Stack[depth].level) &&
                  (Stack[depth - 3].level >= Stack[depth - 1].level) &&
                  (Stack[depth - 3].level >= Stack[depth - 2].level)) {
                if (ok && (!opFlag) && (Stack[depth - 3].type == STYP_SEL3) &&
                    (Stack[depth].type == STYP_VALU) &&
                    (Stack[depth - 1].type == STYP_VALU) &&
                    (Stack[depth - 2].type == STYP_VALU)) {
                  /* 2 argument logical operator */
                  p_return_if_error(
                      SelectorSelect3(G, &Stack[depth - 3], state));
                  opFlag = true;
                  for (a = depth + 1; a <= totDepth; a++)
                    Stack[a - 3] = std::move(Stack[a]);
                  totDepth -= 3;
                }
              }
          if (ok)
            if (depth > 4)
              if ((Stack[depth - 3].level >= Stack[depth].level) &&
                  (Stack[depth - 3].level >= Stack[depth - 1].level) &&
                  (Stack[depth - 3].level >= Stack[depth - 2].level) &&
                  (Stack[depth - 3].level >= Stack[depth - 4].level)) {
                if (ok && (!opFlag) && (Stack[depth - 3].type == STYP_OP22) &&
                    (Stack[depth - 1].type == STYP_VALU) &&
                    (Stack[depth - 2].type == STYP_VALU) &&
                    (Stack[depth].type == STYP_LIST) &&
                    (Stack[depth - 4].type == STYP_LIST)) {
                  ok = SelectorOperator22(G, &Stack[depth - 4], state);
                  opFlag = true;
                  for (a = depth + 1; a <= totDepth; a++)
                    Stack[a - 4] = std::move(Stack[a]);
                  totDepth -= 4;
                }
              }
          if (opFlag) {
            depth = 1; /* start back at the left hand side */
          } else {
            depth = depth + 1;
            opFlag = true;
            if (depth > totDepth) break;
          }
        }
        if (level)
          level--;
        else
          break;
      }
    depth = totDepth;
  }

  if (!ok) {
    return pymol::make_error(indicate_last_token(word, c));
  }

  if (depth != 1) {
    return pymol::make_error("Malformed selection.");
  }

  if (Stack[depth].type != STYP_LIST) {
    return pymol::make_error("Invalid selection.");
  }

  return std::move(Stack[totDepth].sele); /* return the selection list */
}
/**
 * Break a selection down into tokens and return them in a vector.
 * E.g. "(name CA+CB)" -> {"(", "name", "CA+CB", ")"}.
 * @param s selection expression to parse
 * @return tokens
 */
std::vector<std::string> SelectorParse(PyMOLGlobals* G, const char* s) {
  int w_flag = false;
  int quote_flag = false;
  char quote_char = '"';
  const char* p = s;
  std::string* q = nullptr;
  std::vector<std::string> r;
  while (*p) {
    if (w_flag) { /* currently in a word, thus q is a valid pointer */
      if (quote_flag) {
        if (*p != quote_char) {
          *q += *p;
        } else {
          quote_flag = false;
          *q += *p;
        }
      } else
        switch (*p) {
          case ' ':
            w_flag = false;
            break;
          case ';': /* special word terminator */
            *q += *p;
            w_flag = false;
            break;
          case '!': /* single words */
          case '&':
          case '|':
          case '(':
          case ')':
          case '%':
            r.emplace_back(1, *p); /* add new word */
            q = &r.back();
            w_flag = false;
            break;
          case '>':
          case '<':
          case '=':
            if (*(p + 1) == '=') { /* handle >=, <=, == */
              r.emplace_back(p, 2);
              p++;
            } else {
              r.emplace_back(1, *p);
            }
            q = &r.back();
            w_flag = false;
            break;
          case '"':
            quote_flag = true;
            *q += *p;
            break;
          default:
            *q += *p;
            break;
        }
    } else { /*outside a word -- q is undefined */

      switch (*p) {
        case '!': /* single words */
        case '&':
        case '|':
        case '(':
        case ')':
        case '%':
          r.emplace_back(1, *p); /* add new word */
          q = &r.back();
          break;
        case '>':
        case '<':
        case '=':
          if (*(p + 1) == '=') { /* handle >=, <=, == */
            r.emplace_back(p, 2);
            p++;
          } else {
            r.emplace_back(1, *p);
          }
          q = &r.back();
          w_flag = false;
          break;
        case ' ':
          break;
        case '"':
          quote_flag = true;
          quote_char = *p;
          w_flag = true;
          r.emplace_back(1, *p); /* add new word */
          q = &r.back();
          break;
        default:
          w_flag = true;
          r.emplace_back(1, *p); /* add new word */
          q = &r.back();
          break;
      }
    }
    p++;
  }

  if (Feedback(G, FB_Selector, FB_Debugging)) {
    for (auto& word : r) {
      fprintf(stderr, "word: %s\n", word.c_str());
    }
  }
  return (r);
}
