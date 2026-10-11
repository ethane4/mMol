/*
 * CIF tokenizer
 *
 * All keys are canonicalized to lowercase
 *
 * (c) 2014 Schrodinger, Inc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <variant>

#include "pymol/io/cif_file.h"
#include "pymol/base/file.h"
#include "pymol/platform/memory_debug.h"
#include <strings.h>

using namespace pymol::cif;

namespace pymol {
namespace _cif_detail {

template <>
const char* raw_to_typed(const char* s) {
  return s;
}
template <>
std::string raw_to_typed(const char* s) {
  return s;
}
template <>
char raw_to_typed(const char* s) {
  return s[0];
}
template <>
int raw_to_typed(const char* s) {
  return atoi(s);
}

/**
 * Convert to floating point number, ignores uncertainty notation
 * 1.23(45)e2 -> 1.23e2
 */
template <>
double raw_to_typed(const char* s) {
  const char *close, *open = strchr(s, '(');
  if (open && (close = strchr(open, ')'))) {
    return atof(std::string(s, open - s).append(close + 1).c_str());
  }
  return atof(s);
}

template <>
float raw_to_typed(const char* s) {
  return static_cast<float>(raw_to_typed<double>(s));
}

}  // namespace _cif_detail

// basic IO and string handling

// Return true if "c" is whitespace or null
static bool iswhitespace0(char c) {
  return strchr(" \t\r\n", c) ? true : false;
}

// Return true if "c" is whitespace
static bool iswhitespace(char c) { return (c && iswhitespace0(c)); }

// Return true if "c" is line feed or carriage return
static bool islinefeed(char c) { return (c == '\r' || c == '\n'); }

// Return true if "c" is line feed or carriage return or null
static bool islinefeed0(char c) { return (!c || islinefeed(c)); }

// Return true if "c" is double or single quote
static bool isquote(char c) { return (c == '"' || c == '\''); }

// FreeBSD name conflict
#ifdef isspecial
#undef isspecial
#endif

// Return true if token is a STAR keyword
static bool isspecial(const char* token) {
  return (token[0] == '_' || strncasecmp("data_", token, 5) == 0 ||
          strncasecmp("save_", token, 5) == 0 ||
          strcasecmp("loop_", token) == 0 || strcasecmp("stop_", token) == 0 ||
          strcasecmp("global_", token) == 0);
}

// convert all chars to lowercase
static void tolowerinplace(char* p) {
  for (; *p; p++) {
    if (*p <= 'Z' && *p >= 'A') *p -= 'Z' - 'z';
  }
}

// CIF stuff

static const cif_array EMPTY_ARRAY(nullptr);

/*
 * Class to store CIF loops. Only for parsing, do not use in any higher level
 * reading functions.
 */
class cif_loop {
 public:
  int ncols;
  int nrows;
  const char** values;

  // methods
  const char* get_value_raw(int row, int col) const;
};

// get table value, return nullptr if indices out of bounds
const char* cif_loop::get_value_raw(int row, int col) const {
  if (row >= nrows) return nullptr;
  return values[row * ncols + col];
}

// get the number of elements in this array
unsigned cif_array::size() const {
  if (auto arr = std::get_if<cif_detail::cif_str_array>(&m_array)) {
    return (arr->col == cif_detail::cif_str_array::NOT_IN_LOOP)
               ? 1
               : arr->pointer.loop->nrows;
  } else if (auto arr = std::get_if<cif_detail::bcif_array>(&m_array)) {
    return arr->m_arr.size();
  }
  return 0;
}

/// Get array value, return nullptr if `pos >= size()` or value in ['.', '?']
const char* cif_detail::cif_str_array::get_value_raw(unsigned pos) const {
  if (col == NOT_IN_LOOP) return (pos > 0) ? nullptr : pointer.value;
  return pointer.loop->get_value_raw(pos, col);
}

// true if all values in ['.', '?']
bool cif_array::is_missing_all() const {
  for (unsigned i = 0, n = size(); i != n; ++i) {
    if (!is_missing(i)) return false;
  }

  return true;
}

/**
 * Get a pointer to array or nullptr if not found
 *
 * Can lookup different aliases, the first one found is returned.
 * Also supports an alias shortcut for the trivial case where mmCIF uses
 * a colon and CIF uses an underscore: `get_arr("_foo?bar")` is identical to
 * `get_arr("_foo.bar", "_foo_bar")`
 *
 * @param key data name, must be lower case
 */
const cif_array* cif_data::get_arr(const char* key) const {
  if (auto data = std::get_if<pymol::cif_detail::cif_str_data>(&m_data)) {
    const auto& dict = data->m_dict;
    const char* p = strchr(key, '?');
    std::remove_reference_t<decltype(dict)>::const_iterator it;

#ifndef NDEBUG
    for (const char* q = key; *q; ++q) {
      assert("key must be lower case" && !('Z' >= *q && *q >= 'A'));
    }
#endif

    // support alias shortcut: '?' matches '.' and '_'
    if (p != nullptr) {
      std::string tmp(key);
      // replace '?' by '.' or '_'
      tmp[p - key] = '.';
      if ((it = dict.find(tmp.c_str())) != dict.end()) return &it->second;
      tmp[p - key] = '_';
      if ((it = dict.find(tmp.c_str())) != dict.end()) return &it->second;
    } else {
      if ((it = dict.find(key)) != dict.end()) return &it->second;
    }
  } else if (auto data = std::get_if<pymol::cif_detail::bcif_data>(&m_data)) {
    const auto& dict = data->m_dict;

    std::string_view keyView(key);
    auto split_key = [](const char c) {
      return c == '.' /*|| c == '_'*/ || c == '?';
    };
    auto splitTokenIt = std::find_if(keyView.begin(), keyView.end(), split_key);
    if (splitTokenIt == keyView.end()) {
      return nullptr;
    }
    auto dist = std::distance(keyView.begin(), splitTokenIt);
    auto categoryView = keyView.substr(0, dist);
    auto categoryStr = std::string(categoryView);
    auto categoryIt = dict.find(categoryStr.c_str());
    if (categoryIt == dict.end()) {
      return nullptr;
    }
    auto& category = categoryIt->second;
    auto columnView = keyView.substr(dist + 1);
    auto columnStr = std::string(columnView);
    auto columnIt = category.find(columnStr.c_str());
    if (columnIt == category.end()) {
      return nullptr;
    }
    return &columnIt->second;
  }

  return nullptr;
}

const char* cif_data::code() const {
  if (auto data = std::get_if<pymol::cif_detail::cif_str_data>(&m_data)) {
    return data->m_code ? data->m_code : "";
  }
  if (auto data = std::get_if<pymol::cif_detail::bcif_data>(&m_data)) {
    return data->m_code.c_str();
  }
  return "";
}

const cif_array* cif_data::empty_array() { return &EMPTY_ARRAY; }

bool cif_file::parse_file(const char* filename) {
  char* contents = FileGetContents(filename, nullptr);

  if (!contents) {
    error(std::string("failed to read file ").append(filename).c_str());
    return false;
  }

  return parse(std::move(contents));
}

bool cif_file::parse_string(const char* contents) {
  return parse(std::move(mstrdup(contents)));
}

void cif_file::error(const char* msg) {
  std::cout << "ERROR " << msg << std::endl;
}

// constructor
cif_file::cif_file(const char* filename, const char* contents_) {
  if (contents_) {
    parse_string(contents_);
  } else if (filename) {
    parse_file(filename);
  }
}

// constructor
cif_file::cif_file() = default;
cif_file::cif_file(cif_file&&) = default;

// move assignment
cif_file& cif_file::operator=(cif_file&&) = default;

// destructor
cif_file::~cif_file() = default;

bool cif_file::parse(char*&& p) {
  m_datablocks.clear();
  m_tokens.clear();
  m_contents.reset(p);

  if (!p) {
    error("parse(nullptr)");
    return false;
  }

  auto& tokens = m_tokens;
  char quote;
  char prev = '\0';

  std::vector<bool> keypossible;

  // tokenize
  while (true) {
    while (iswhitespace(*p)) prev = *(p++);

    if (!*p) break;

    if (*p == '#') {
      while (!(islinefeed0(*++p)));
      prev = *p;
    } else if (isquote(*p)) {  // will nullptr the closing quote
      quote = *p;
      keypossible.push_back(false);
      tokens.push_back(p + 1);
      while (*++p && !(*p == quote && iswhitespace0(p[1])));
      if (*p) *(p++) = 0;
      prev = *p;
    } else if (*p == ';' && islinefeed(prev)) {
      // multi-line tokens start with ";" and end with "\n;"
      // multi-line tokens cannot be keys, only values.
      keypossible.push_back(false);
      tokens.push_back(p + 1);
      // advance until `\n;`
      while (*++p && !(islinefeed(*p) && p[1] == ';'));
      // step to next line and null the line feed
      if (*p) {
        *p = 0;
        // \r\n on Windows)
        if (p - 1 > tokens.back() && *(p - 1) == '\r') {
          *(p - 1) = 0;
        }
        p += 2;
      }
      prev = ';';
    } else {  // will null the whitespace
      char* q = p++;
      while (!iswhitespace0(*p)) ++p;
      prev = *p;
      if (p - q == 1 && (*q == '?' || *q == '.')) {
        // store values '.' (inapplicable) and '?' (unknown) as null-pointers
        q = nullptr;
        keypossible.push_back(false);
      } else {
        if (*p) *(p++) = 0;
        keypossible.push_back(true);
      }
      tokens.push_back(q);
    }
  }

  cif_detail::cif_str_data* current_frame = nullptr;
  std::vector<cif_detail::cif_str_data*> frame_stack;
  std::unique_ptr<cif_data> global_block;
  decltype(m_datablocks) datablocksnew;

  // parse into dictionary
  for (unsigned int i = 0, n = tokens.size(); i < n; i++) {
    if (!keypossible[i]) {
      error("expected key (1)");
      return false;
    } else if (tokens[i][0] == '_') {
      if (!current_frame) {
        error("missing data_ (unexpected data name)");
        return false;
      }

      if (i + 1 == n) {
        error("truncated");
        return false;
      }

      tolowerinplace(tokens[i]);
      current_frame->m_dict[tokens[i]].m_array = cif_detail::cif_str_array{};
      auto& cif_arr = std::get<cif_detail::cif_str_array>(
          current_frame->m_dict[tokens[i]].m_array);
      cif_arr.set_value(tokens[i + 1]);

      i++;
    } else if (strcasecmp("loop_", tokens[i]) == 0) {
      if (!current_frame) {
        error("missing data_ (unexpected loop)");
        return false;
      }

      int ncols = 0;
      int nrows = 0;
      cif_loop* loop = nullptr;

      // loop data
      loop = new cif_loop;
      current_frame->m_loops.emplace_back(loop);

      // columns
      while (++i < n && keypossible[i] && tokens[i][0] == '_') {
        tolowerinplace(tokens[i]);
        current_frame->m_dict[tokens[i]].m_array = cif_detail::cif_str_array{};
        auto& cif_arr = std::get<cif_detail::cif_str_array>(
            current_frame->m_dict[tokens[i]].m_array);
        cif_arr.set_loop(loop, ncols);

        ncols++;
      }

      if (loop) {
        // loop data
        loop->values = (const char**)&tokens[i];
        loop->ncols = ncols;
      }

      // rows
      while (i < n && !(keypossible[i] && isspecial(tokens[i]))) {
        i += ncols;

        if (i > n) {
          error("truncated loop");
          return false;
        }

        nrows++;
      }

      // loop data
      if (loop) {
        loop->nrows = nrows;
      }

      i--;

    } else if (strncasecmp("data_", tokens[i], 5) == 0) {
      auto& new_data = datablocksnew[tokens[i] + 5];
      new_data.m_data = cif_detail::cif_str_data();
      current_frame = &std::get<cif_detail::cif_str_data>(new_data.m_data);
      current_frame->m_code = tokens[i] + 5;
      frame_stack = {current_frame};

    } else if (strncasecmp("global_", tokens[i], 5) == 0) {
      // STAR feature, not supported in CIF
      auto new_data = new cif_data;
      new_data->m_data = cif_detail::cif_str_data{};
      current_frame = &std::get<cif_detail::cif_str_data>(new_data->m_data);
      global_block.reset(new_data);
      frame_stack = {current_frame};

    } else if (strncasecmp("save_", tokens[i], 5) == 0) {
      if (tokens[i][5]) {
        // begin
        if (!current_frame) {
          error("top-level save_");
          return false;
        }

        const char* key(tokens[i] + 5);
        current_frame = &current_frame->m_saveframes[key];
        frame_stack.push_back(current_frame);
      } else {
        // end
        if (frame_stack.size() < 2) {
          error("unexpected save_");
          return false;
        }

        frame_stack.pop_back();
        current_frame = frame_stack.back();
      }
    } else {
      error("expected key (2)");
      return false;
    }
  }

  m_datablocks = std::move(datablocksnew);

  return true;
}

bool cif_file::parse_bcif(const char* bytes, std::size_t size) { return false; }

}  // namespace pymol
