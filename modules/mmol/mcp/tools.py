"""MCP tool definitions + dispatcher for the mMol MCP server.

Tools call PyMOL's ``cmd`` API, imported LAZILY inside functions so this module
imports standalone (for unit tests) without a built PyMOL. Handlers run on the
MCP server's HTTP request-handler threads, but every ``cmd``-touching body is
run via ``mainthread.run_on_main`` so it executes on the app's MAIN thread — the
only thread where the embedded core is safe to touch (it races
``SceneRenderMetal`` otherwise; see mainthread.py). Never call the C bridge
(PyMOLBridge_*) from here; only the Python ``cmd`` API.

``capture_viewport`` renders with ``ray`` (cmd.png ray=1; GPU when it can)
because cmd.png without ray reads a GL framebuffer the Metal app does not have.
"""

import base64
import contextlib
import io
import json
import os
import tempfile
import threading
import traceback
import urllib.parse
import urllib.request

from mmol.mcp.mainthread import run_on_main

# Persistent namespace for run_python; seeded on first use.
_py_ns = None
_py_ns_lock = threading.Lock()


def _namespace():
    global _py_ns
    with _py_ns_lock:
        if _py_ns is None:
            from pymol import cmd
            ns = {"__name__": "__mcp_exec__", "cmd": cmd}
            try:
                import numpy as np
                ns["np"] = np
            except Exception:
                pass
            try:
                import Bio
                ns["Bio"] = Bio
            except Exception:
                pass
            _py_ns = ns
    return _py_ns


def _text(s):
    return {"content": [{"type": "text", "text": str(s)}], "isError": False}


def _error(s):
    return {"content": [{"type": "text", "text": str(s)}], "isError": True}


# --- Tool implementations -------------------------------------------------


def _run_pymol_command(args):
    cmd_str = args.get("command", "")
    if not cmd_str:
        return _error("missing 'command'")

    def _body():
        from pymol import cmd
        cmd.do(cmd_str)
        cmd.sync()
        return _text("ok")

    return run_on_main(_body)


def _run_python(args):
    code = args.get("code", "")
    if not code:
        return _error("missing 'code'")

    def _body():
        ns = _namespace()
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                exec(code, ns)
            from pymol import cmd
            cmd.sync()
            out = buf.getvalue()
            return _text(out if out else "ok")
        except Exception:
            return _error(buf.getvalue() + "\n" + traceback.format_exc())

    return run_on_main(_body)


def _get_session_state(args):

    def _body():
        from pymol import cmd
        objects = []
        for name in cmd.get_names("objects"):
            otype = cmd.get_type(name)
            # Only molecule objects are valid atom-selections. A measurement,
            # map, CGO or group name fed to count_atoms('(name)') raises a C++
            # "Invalid selection name" Selector-Error (same gotcha mmol/pick.py
            # documents), which would otherwise abort the whole call.
            objects.append({
                "name":
                    name,
                "type":
                    otype,
                "atoms":
                    (cmd.count_atoms("(%s)" %
                                     name) if otype == "object:molecule" else 0
                    ),
            })
        state = {
            "objects": objects,
            "selections": cmd.get_names("selections"),
            "view": list(cmd.get_view()),
            "frame": cmd.get_frame(),
            "n_frames": cmd.count_frames(),
        }
        return _text(json.dumps(state, indent=2))

    return run_on_main(_body)


def _capture_viewport(args):
    path = None
    try:
        # Clamp to a sane maximum: width/height come straight from the caller and
        # feed a ray-trace, so an unbounded size is a multi-minute / OOM DoS.
        width = max(1, min(int(args.get("width", 640)), 2048))
        height = max(1, min(int(args.get("height", 480)), 2048))
        fd, path = tempfile.mkstemp(suffix=".png", prefix="mmol_mcp_capture_")
        os.close(fd)

        # ray=1 ray-traces and writes the PNG in one synchronous call. A separate
        # cmd.ray + cmd.png(prior=1) fails with "no prior image available" when
        # driven from the MCP server thread. Runs on the main thread (touches the
        # scene/renderer); the tempfile create/read/cleanup stay off-main.
        def _render():
            from pymol import cmd
            cmd.png(path, width=width, height=height, ray=1)

        run_on_main(_render)

        with open(path, "rb") as f:
            data = base64.b64encode(f.read()).decode("ascii")
        return {
            "content": [{
                "type": "image",
                "data": data,
                "mimeType": "image/png"
            }],
            "isError": False
        }
    except Exception:
        return _error(traceback.format_exc())
    finally:
        if path:
            try:
                os.remove(path)
            except OSError:
                pass


def _search_pdb(args):
    query = args.get("query", "")
    if not query:
        return _error("missing 'query'")
    try:
        limit = int(args.get("limit", 10))
        body = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": query
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": limit
                }
            },
        }
        url = ("https://search.rcsb.org/rcsbsearch/v2/query?json=" +
               urllib.parse.quote(json.dumps(body)))
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        ids = [r["identifier"] for r in data.get("result_set", [])]
        return _text(json.dumps({"pdb_ids": ids}, indent=2))
    except Exception:
        return _error(traceback.format_exc())


# --- mMol additions: CSV heatmaps and Cα contact maps ---------------------


def _emit(kind, payload):
    """Hand a request to the Swift app on the feedback line (polled every
    100 ms), base64 JSON so it survives the line split. Main thread only."""
    print("%s:%s" %
          (kind, base64.b64encode(json.dumps(payload).encode()).decode()))


def _color_by_csv(args):
    obj = args.get("object", "")
    path = args.get("csv_path") or ""
    text = args.get("csv_text") or ""
    if not obj or not (path or text):
        return _error("need 'object' and one of 'csv_path' / 'csv_text'")
    if text and not path:
        fd, path = tempfile.mkstemp(prefix="mmol_heatmap_", suffix=".csv")
        with os.fdopen(fd, "w") as f:
            f.write(text)
    path = os.path.expanduser(path)
    try:
        with open(path) as f:
            header = f.readline()
            nrows = sum(1 for line in f if line.strip())
    except OSError as e:
        return _error("can't read %s: %s" % (path, e))
    sep = "," if ("," in header or "\t" not in header) else "\t"
    cols = [c.strip().strip('"') for c in header.strip().split(sep)]
    for key in ("column", "index_column"):
        v = args.get(key)
        if v and v not in cols:
            return _error("no column %r in the CSV; columns are: %s" %
                          (v, ", ".join(cols)))

    def _rgb(v):
        if v is None:
            return None
        if isinstance(v, (list, tuple)) and len(v) == 3:
            return [float(x) for x in v]
        s = str(v).strip()
        if s.startswith("#") and len(s) == 7:
            return [int(s[i:i + 2], 16) / 255.0 for i in (1, 3, 5)]
        from pymol import cmd
        return list(cmd.get_color_tuple(s))

    def _body():
        from pymol import cmd
        if obj not in cmd.get_names("objects"):
            return _error("no object named %r" % obj)
        req = {"object": obj, "path": path}
        for key in ("column", "index_column", "min", "max", "through_white"):
            if key in args:
                req[key] = args[key]
        for key in ("low_color", "high_color"):
            if args.get(key) is not None:
                req[key] = _rgb(args[key])
        _emit("HEATMAP", req)
        return _text(
            "Heatmap requested for %s from %s (%d rows; columns: %s). The app "
            "colours the object and shows the settings in its \u25b8 card's "
            "heatmap chip; residues the CSV doesn't cover turn grey." %
            (obj, os.path.basename(path), nrows, ", ".join(cols)))

    return run_on_main(_body)


def _contact_map(args):
    obj = args.get("object", "")
    if not obj:
        return _error("missing 'object'")
    cutoff = float(args.get("cutoff", 5.0))
    min_sep = int(args.get("min_seq_separation", 1))
    max_pairs = int(args.get("max_pairs", 500))
    open_tab = bool(args.get("open_tab", True))

    def _body():
        from pymol import cmd
        if obj not in cmd.get_names("objects"):
            return _error("no object named %r" % obj)
        sel = "(%s) and polymer and name CA and present and (alt '' or alt A)" % obj
        res, xyz = [], []
        cmd.iterate_state(
            -1,
            sel,
            "res.append((chain, resi, resn)); xyz.append((x, y, z))",
            space={
                "res": res,
                "xyz": xyz
            })
        if not res:
            return _error("%s has no protein C-alpha atoms" % obj)
        label = lambda k: "%s%s%s" % (res[k][0] + ":" if res[k][0] else "", res[
            k][2], res[k][1])
        sep = max(1, min_sep)
        try:
            import numpy as np
            c = np.asarray(xyz, dtype=float)
            d = np.sqrt(((c[:, None, :] - c[None, :, :])**2).sum(-1))
            i, j = np.nonzero(np.triu(d <= cutoff, k=sep))
            pairs = [[label(a), label(b),
                      round(float(d[a, b]), 2)] for a, b in zip(i, j)]
        except ImportError:
            # Builds without numpy: the same pairs, in plain Python.
            c2, n, pairs = cutoff * cutoff, len(xyz), []
            for a in range(n):
                xa, ya, za = xyz[a]
                for b in range(a + sep, n):
                    xb, yb, zb = xyz[b]
                    q = (xa - xb)**2 + (ya - yb)**2 + (za - zb)**2
                    if q <= c2:
                        pairs.append([label(a), label(b), round(q**0.5, 2)])
        if open_tab:
            _emit("CONTACTMAP", {"object": obj})
        return _text(
            json.dumps(
                {
                    "object": obj,
                    "cutoff_A": cutoff,
                    "residues": len(res),
                    "min_seq_separation": min_sep,
                    "contacts": len(pairs),
                    "pairs": pairs[:max_pairs],
                    "truncated": len(pairs) > max_pairs,
                    "tab_opened": open_tab
                },
                indent=1))

    return run_on_main(_body)


# --- Registry + dispatcher ------------------------------------------------

TOOLS = [
    {
        "name": "run_pymol_command",
        "description": ("Run one PyMOL command-language statement (e.g. "
                        "'fetch 1ubq, async=0', 'show cartoon', "
                        "'color red, chain A'). For multi-step logic or data "
                        "access, prefer run_python."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "command": {
                    "type": "string"
                }
            },
            "required": ["command"],
        },
    },
    {
        "name": "run_python",
        "description":
            ("Execute arbitrary Python in a persistent namespace with "
             "'cmd' (the PyMOL API), and 'np'/'Bio' if available. "
             "State persists across calls. stdout is returned. This is "
             "the most powerful tool -- prefer it for anything beyond a "
             "single command."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "code": {
                    "type": "string"
                }
            },
            "required": ["code"],
        },
    },
    {
        "name": "get_session_state",
        "description": ("Return the current session as JSON: objects "
                        "(name/type/atom count), named selections, camera view "
                        "(18 floats), current frame and frame count."),
        "inputSchema": {
            "type": "object",
            "properties": {}
        },
    },
    {
        "name": "capture_viewport",
        "description":
            ("Ray-trace the current view (PyMOL's `ray`, on the GPU "
             "when it can) and return it as a PNG image so you can see "
             "the structure. Optional width/height (default 640x480)."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "width": {
                    "type": "integer",
                    "default": 640
                },
                "height": {
                    "type": "integer",
                    "default": 480
                },
            },
        },
    },
    {
        "name": "search_pdb",
        "description":
            ("Full-text search of the RCSB PDB; returns matching PDB "
             "IDs. Use run_pymol_command 'fetch <id>, async=0' to load "
             "one."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string"
                },
                "limit": {
                    "type": "integer",
                    "default": 10
                },
            },
            "required": ["query"],
        },
    },
    {
        "name": "color_by_csv",
        "description": (
            "Colour an object as a heatmap from a CSV column (mMol's "
            "C \u25b8 by file\u2026). The CSV has a header row and one row per "
            "residue: with no index column, row 1 is residue 1; an index "
            "column's values are residue numbers ('45', every chain) or "
            "chain-qualified ('A:45' / 'A45'). Pass csv_path (a file on this "
            "Mac) or csv_text (the CSV itself). Optional: column (default "
            "the first numeric column), index_column ('' for none; a first "
            "column named resi/residue/index is picked automatically), "
            "low_color / high_color (PyMOL name, '#rrggbb' or [r,g,b] 0-1; "
            "default blue/red), through_white (default true: low \u2192 white "
            "\u2192 high), min / max (default the column's range). Uncovered "
            "residues turn grey. The settings appear in the object's "
            "heatmap panel, where the user can keep adjusting them."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "object": {
                    "type": "string"
                },
                "csv_path": {
                    "type": "string"
                },
                "csv_text": {
                    "type": "string"
                },
                "column": {
                    "type": "string"
                },
                "index_column": {
                    "type": "string"
                },
                "low_color": {},
                "high_color": {},
                "through_white": {
                    "type": "boolean"
                },
                "min": {
                    "type": "number"
                },
                "max": {
                    "type": "number"
                },
            },
            "required": ["object"],
        },
    },
    {
        "name": "contact_map",
        "description":
            ("C-alpha contact map of an object (mMol's A \u25b8 Contact "
             "Map\u2026): residue pairs whose C-alpha atoms are within "
             "cutoff \u00c5 (default 5) in the current state, as JSON "
             "[residue, residue, distance]. min_seq_separation (default 1) "
             "drops pairs closer in sequence than that (e.g. 3 to skip "
             "i,i+1 and i,i+2 neighbours). max_pairs caps the list (default "
             "500; the total count is always given). open_tab (default "
             "true) also opens the interactive map tab for the user."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "object": {
                    "type": "string"
                },
                "cutoff": {
                    "type": "number",
                    "default": 5
                },
                "min_seq_separation": {
                    "type": "integer",
                    "default": 1
                },
                "max_pairs": {
                    "type": "integer",
                    "default": 500
                },
                "open_tab": {
                    "type": "boolean",
                    "default": True
                },
            },
            "required": ["object"],
        },
    },
]

_DISPATCH = {
    "color_by_csv": _color_by_csv,
    "contact_map": _contact_map,
    "run_pymol_command": _run_pymol_command,
    "run_python": _run_python,
    "get_session_state": _get_session_state,
    "capture_viewport": _capture_viewport,
    "search_pdb": _search_pdb,
}


def call(name, arguments):
    fn = _DISPATCH.get(name)
    if fn is None:
        return _error("unknown tool: %s" % name)
    try:
        return fn(arguments or {})
    except Exception:
        return _error(traceback.format_exc())
