"""`ray`: ray-traced images shown in the viewport like PyMOL 2.5.

Active only in the Classic PyMOL Layout (the Swift side calls enable()).

  ray   [width, height]   PyMOL's own ray tracer (run on the GPU when it can, same
                          image), as in desktop PyMOL 2.5 — every
                          ray_* / lighting setting (e.g. from ~/.pymolrc) applies:
                          ray_trace_mode outlines, ray_trace_color, ray_shadows,
                          ray_trace_gain, specular, ambient, direct, light …

The image ends up where PyMOL keeps a ray result (the scene image), so
a following `png file.png` saves it and any change to the scene drops it; and a
RAYIMG: feedback marker makes the viewport show it over the live scene
(checkerboard behind transparency, outlined when smaller than the viewport) until
the core drops it. The Metal viewport cannot draw PyMOL's own image overlay.
"""

import os
import struct
import tempfile

from pymol import cmd

_enabled = False
_orig_ray = None


def _bg_rgb(_self):
    """bg_rgb as floats. cmd.get gives a colour name ("black") or, for an RGB
    value, a "[ r, g, b ]" string."""
    val = _self.get('bg_rgb')
    try:
        if val.strip().startswith('['):
            return tuple(float(x) for x in val.strip('[] ').split(','))[:3]
        return tuple(_self.get_color_tuple(val))[:3]
    except Exception:
        return (0.0, 0.0, 0.0)


def _png_size(path):
    with open(path, 'rb') as f:
        head = f.read(24)
    return struct.unpack('>II', head[16:24])


def _show(_self, path):
    """Tell the viewport to show the scene image saved at `path`."""
    w, h = _png_size(path)
    r, g, b = _bg_rgb(_self)
    print("RAYIMG:%d:%d:%.3f,%.3f,%.3f:%s" % (w, h, r, g, b, path))


def _tmp_png():
    fd, path = tempfile.mkstemp(prefix='mmol_ray_', suffix='.png')
    os.close(fd)
    return path


def _ray(width=0,
         height=0,
         antialias=-1,
         angle=0.0,
         shift=0.0,
         renderer=-1,
         quiet=1,
         async_=0,
         _self=cmd,
         **kwargs):
    """PyMOL's ray tracer, then show the result in the Metal viewport."""
    async_ = int(kwargs.pop('async', async_))
    r = _orig_ray(width,
                  height,
                  antialias,
                  angle,
                  shift,
                  renderer,
                  quiet,
                  0 if _enabled else async_,
                  _self=_self,
                  **kwargs)
    if _enabled and not float(angle) and not float(shift):
        import _mmol_render
        path = _tmp_png()
        # a throwaway file the viewport reads back: skip compression
        _mmol_render.fast_png(1)
        try:
            # prior=1: write the image the ray tracer just made, no re-render.
            _self.png(path, prior=1, quiet=1)
            _show(_self, path)
        except Exception as e:
            print(" Ray: couldn't show the image (%s)." % (e,))
        finally:
            _mmol_render.fast_png(0)
    return r


# --- Pale default colours ---------------------------------------------------
#
# PyMOL gives each newly loaded object's carbons the next colour of its auto
# cycle (green, cyan, lightmagenta, yellow, …). In the Classic PyMOL Layout a
# freshly loaded structure gets a pale tint of that colour instead (green →
# pale green): the colour mixed _TINT of the way to white, defined as
# "tint_<name>"; its organic molecules (ligands, cofactors) get a pale orange
# instead. Only an object whose carbons all carry one colour — i.e. still the
# auto colour — is touched, so sessions and files with their own colours keep
# them.

_TINT = 0.55
_orig_load = None
_orig_fetch = None


def _pale(idx, _self=cmd):
    """The load tint of color `idx`: defined once as tint_<idx>, returned."""
    tname = 'tint_%s' % idx
    rgb = _self.get_color_tuple(idx)
    _self.set_color(tname, [c + (1.0 - c) * _TINT for c in rgb], quiet=1)
    return tname


_orig_cbc = None


def _cbc(selection='(all)', first_color=7, quiet=1, legacy=0, _self=cmd):
    """util.cbc (C > by chain) in the same pale tints as new structures."""
    from pymol import util
    if not _enabled or int(legacy):
        return _orig_cbc(selection, first_color, quiet, legacy, _self=_self)
    for c, a in enumerate(_self.get_chains(selection)):
        if len(a.split()) != 1:
            a = '"' + a + '"'
        color = _pale(util._color_cycle[c % util._color_cycle_len], _self)
        _self.color(color, "(chain %s and (%s))" % (a, selection), quiet=quiet)


def _tint_new(before, _self=cmd):
    for obj in _self.get_names('objects'):
        if obj in before or _self.get_type(obj) != 'object:molecule':
            continue
        cols = set()
        _self.iterate('(%s) and elem C' % obj,
                      'cols.add(color)',
                      space={'cols': cols})
        if len(cols) != 1:
            continue
        idx = cols.pop()
        try:
            rgb = _self.get_color_tuple(idx)
            name = str(idx)
        except Exception:
            continue
        tint = [c + (1.0 - c) * _TINT for c in rgb]
        tname = 'tint_%s' % name
        _self.set_color(tname, tint, quiet=1)
        _self.color(tname, '(%s) and elem C' % obj, quiet=1)
        # ligands / cofactors: one fixed contrasting tint (pale orange)
        _self.set_color('tint_ligand',
                        [c + (1.0 - c) * _TINT for c in (1.0, 0.5, 0.0)],
                        quiet=1)
        _self.color('tint_ligand', '(%s) and organic and elem C' % obj, quiet=1)


def _wrap_loader(orig):
    import functools

    @functools.wraps(
        orig)  # keeps the doc and signature (`help`, console hints)
    def loader(*args, **kwargs):
        if not _enabled:
            return orig(*args, **kwargs)
        target = str(
            args[0] if args else kwargs.get('filename', kwargs.get('code', '')))
        if orig is _orig_fetch:
            # Wait for the download, so the new object exists to be tinted.
            kwargs.pop('async', None)
            kwargs['async_'] = 0
        before = set(cmd.get_names('objects'))
        r = orig(*args, **kwargs)
        if not target.lower().endswith(('.pse', '.psw')):
            try:
                _tint_new(before)
            except Exception:
                pass
        return r

    return loader


def enable(on=True):
    """Classic layout on/off: `ray` shows its image in the viewport and new
    structures load in pale tints."""
    global _enabled, _orig_ray, _orig_load, _orig_fetch, _orig_cbc
    if _orig_ray is None:
        from pymol import util
        _orig_cbc = util.cbc
        _cbc.__doc__ = _orig_cbc.__doc__
        util.cbc = _cbc
        if 'util.cbc' in cmd.keyword:
            cmd.keyword['util.cbc'][0] = _cbc
        _orig_ray = cmd.ray
        cmd.ray = _ray
        _ray.__doc__ = _orig_ray.__doc__  # `help ray`, console hints
        cmd.keyword['ray'][0] = _ray
        _orig_load, _orig_fetch = cmd.load, cmd.fetch
        cmd.load = _wrap_loader(_orig_load)
        cmd.keyword['load'][0] = cmd.load
        cmd.fetch = _wrap_loader(_orig_fetch)
        cmd.keyword['fetch'][0] = cmd.fetch
    _enabled = bool(on)
