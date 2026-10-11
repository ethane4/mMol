"""Regression test for mMol, run inside the app by ``run.sh``.

Renders a few scenes with the GPU and CPU ray tracers and records the default
settings and the results of core commands (symmetry, alignment, measurement,
scenes, movies, sessions) and of the app's ``mmol`` modules. Everything is
compared with the reference files in ``golden/``; with ``MMOL_TEST_UPDATE=1``
the references are rewritten instead.

The structures load from ``data/``, and saved camera views from
``golden/views.json``, so the results do not depend on the network or on the
window size.

Environment
-----------
MMOL_TEST_ROOT
    This directory.
MMOL_TEST_OUT
    Where to write the renders, the session file and ``report.txt``.
MMOL_TEST_UPDATE
    ``1`` to rewrite the reference files.
"""

import importlib
import json
import os
import shutil
import struct
import traceback
import zlib

import _mmol_render
from pymol import cmd

ROOT = os.environ['MMOL_TEST_ROOT']
OUT = os.environ['MMOL_TEST_OUT']
UPDATE = os.environ.get('MMOL_TEST_UPDATE') == '1'
GOLDEN = os.path.join(ROOT, 'golden')
DATA = os.path.join(ROOT, 'data')
VIEWS = {}  # camera views by scene name, from golden/views.json


def rounded(value, digits=4):
    """Round the floats in a JSON-like value so the result is stable.

    Parameters
    ----------
    value : object
        A number, string, or nested lists, tuples and dicts of them.
    digits : int, optional
        Decimal places to keep.

    Returns
    -------
    object
        `value` with every float rounded and tuples turned into lists.
    """
    if isinstance(value, float):
        return round(value, digits)
    if isinstance(value, (list, tuple)):
        return [rounded(v, digits) for v in value]
    if isinstance(value, dict):
        return {k: rounded(v, digits) for k, v in value.items()}
    return value


def read_png(path):
    """Decode an 8-bit grayscale, RGB or RGBA PNG file.

    Parameters
    ----------
    path : str
        The PNG file.

    Returns
    -------
    tuple of (int, int, bytes)
        Width, height and the unfiltered pixel bytes, row after row.
    """
    data = open(path, 'rb').read()
    pos, idat = 8, b''
    while pos < len(data):
        size, kind = struct.unpack('>I4s', data[pos:pos + 8])
        body = data[pos + 8:pos + 8 + size]
        if kind == b'IHDR':
            width, height, depth, color = struct.unpack('>IIBB', body[:10])
        elif kind == b'IDAT':
            idat += body
        pos += 12 + size
    assert depth == 8, 'only 8-bit PNG files are supported'
    bpp = {0: 1, 2: 3, 4: 2, 6: 4}[color]
    stride = width * bpp
    raw = zlib.decompress(idat)
    out = bytearray()
    prev = bytearray(stride)
    for y in range(height):
        start = y * (stride + 1)
        kind, cur = raw[start], bytearray(raw[start + 1:start + 1 + stride])
        for x in range(stride):
            a = cur[x - bpp] if x >= bpp else 0
            b = prev[x]
            c = prev[x - bpp] if x >= bpp else 0
            if kind == 1:
                pred = a
            elif kind == 2:
                pred = b
            elif kind == 3:
                pred = (a + b) // 2
            elif kind == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                pred = a if pa <= pb and pa <= pc else b if pb <= pc else c
            else:
                break
            cur[x] = (cur[x] + pred) & 255
        out += cur
        prev = cur
    return width, height, bytes(out)


class Report:
    """Collects the checks and writes ``report.txt``."""

    def __init__(self):
        self.lines = []
        self.failures = 0

    def add(self, ok, name, detail=''):
        """Record one check.

        Parameters
        ----------
        ok : bool
            Whether the check passed.
        name : str
            What was checked.
        detail : str, optional
            Why it failed.
        """
        self.failures += not ok
        self.lines.append(
            '%s %s%s' %
            ('ok  ' if ok else 'FAIL', name, ': ' + detail if detail else ''))

    def write(self):
        """Write the report, ending with a ``DONE`` line."""
        self.lines.append('DONE %d failure(s)' % self.failures)
        with open(os.path.join(OUT, 'report.txt'), 'w') as f:
            f.write('\n'.join(self.lines) + '\n')


def compare_json(report, name, new):
    """Compare a JSON-able value with ``golden/<name>.json``.

    Parameters
    ----------
    report : Report
        Where to record the result.
    name : str
        The reference file's base name.
    new : dict
        The values recorded in this run.
    """
    path = os.path.join(GOLDEN, name + '.json')
    if UPDATE:
        with open(path, 'w') as f:
            json.dump(new, f, indent=1, sort_keys=True)
            f.write('\n')
        return
    old = json.load(open(path))
    bad = [k for k in sorted(set(old) | set(new)) if old.get(k) != new.get(k)]
    for key in bad:
        report.add(False, '%s %s' % (name, key),
                   'expected %r, got %r' % (old.get(key), new.get(key)))
    if not bad:
        report.add(True, '%s (%d values)' % (name, len(old)))


def compare_image(report, name):
    """Compare ``OUT/<name>`` with ``golden/<name>`` pixel by pixel.

    Parameters
    ----------
    report : Report
        Where to record the result.
    name : str
        The PNG file name.
    """
    new, old = os.path.join(OUT, name), os.path.join(GOLDEN, name)
    if UPDATE:
        shutil.copy(new, old)
        return
    if open(new, 'rb').read() == open(old, 'rb').read():
        report.add(True, name)
        return
    a, b = read_png(old), read_png(new)
    if a[:2] != b[:2]:
        report.add(False, name, 'size %dx%d, expected %dx%d' % (b[:2] + a[:2]))
        return
    diffs = [abs(x - y) for x, y in zip(a[2], b[2]) if x != y]
    report.add(
        not diffs, name, '%.3f%% of values differ, by up to %d' %
        (100.0 * len(diffs) / len(a[2]), max(diffs, default=0)))


def check(results, name, fn, keep=True):
    """Run one step and record its value, or its error.

    Parameters
    ----------
    results : dict
        Maps step names to recorded values.
    name : str
        The step's name.
    fn : callable
        Runs the step and returns a JSON-able value.
    keep : bool, optional
        False to record only that the step ran, for values that depend on the
        window size.
    """
    try:
        value = fn()
        results[name] = rounded(value) if keep else 'ran'
    except Exception:
        results[name] = 'ERROR ' + traceback.format_exc().splitlines()[-1]


def load(name):
    """Load ``data/<name>.cif.gz`` as object `name`."""
    cmd.load(os.path.join(DATA, name + '.cif.gz'), name)


def render(views, name, width, height, setup):
    """Ray trace a scene with the GPU and the CPU ray tracers.

    Parameters
    ----------
    views : dict
        Saved camera views by scene name; updated when recording.
    name : str
        The scene name; the images are ``<name>_gpu.png`` and
        ``<name>_cpu.png``.
    width, height : int
        Image size in pixels.
    setup : callable
        Loads and styles the scene, and orients the camera.
    """
    cmd.delete('all')
    setup()
    if UPDATE:
        views[name] = list(cmd.get_view())
    for kind, gpu in (('gpu', 1), ('cpu', 0)):
        _mmol_render.set_gpu_ray(gpu)
        cmd.set_view(views[name])
        cmd.ray(width, height)
        cmd.png(os.path.join(OUT, '%s_%s.png' % (name, kind)), prior=1, quiet=1)
    _mmol_render.set_gpu_ray(1)


def protease():
    """HIV protease with its inhibitor: cartoon, sticks and a surface."""
    load('1hsg')
    cmd.hide('everything')
    cmd.show('cartoon')
    cmd.show('sticks', 'organic')
    cmd.show('surface', 'chain A and resi 1-30')
    cmd.orient(animate=0)


def ubiquitin():
    """Ubiquitin with spheres, dots, a surface and labels in every font."""
    load('1ubq')
    cmd.hide('everything')
    cmd.show('cartoon')
    cmd.show('spheres', 'resi 1-6')
    cmd.show('dots', 'resi 20-24')
    cmd.show('surface', 'resi 40-50')
    cmd.show('nb_spheres', 'resn HOH')
    for i, font in enumerate(range(5, 19)):
        name = 'font%d' % font
        cmd.create(name, '1ubq and resi %d and name CA' % (55 + i))
        cmd.hide('everything', name)
        cmd.label(name, '"Font %d Ag"' % font)
        cmd.set('label_font_id', font, name)
        cmd.set('label_size', 20, name)
    cmd.orient('1ubq', animate=0)
    cmd.turn('y', 20)


def ramp():
    """Ubiquitin's surface colored by a map ramp, in the CMYK color space."""
    load('1ubq')
    cmd.remove('solvent')
    cmd.hide('everything')
    cmd.map_new('map', 'gaussian', 1.0, '1ubq', 3)
    cmd.ramp_new('ramp', 'map', [0.0, 0.5, 1.0], ['blue', 'white', 'red'])
    cmd.show('surface')
    cmd.set('surface_color', 'ramp', '1ubq')
    cmd.space('cmyk', quiet=1)
    cmd.orient('1ubq', animate=0)


def dna():
    """B-DNA as a cartoon with base rings."""
    load('1bna')
    cmd.remove('solvent')
    cmd.hide('everything')
    cmd.set('cartoon_ring_mode', 3)
    cmd.show('cartoon')
    cmd.orient(animate=0)


def settings():
    """Return every setting's value; run before anything changes them."""
    out = {}
    for name in cmd.setting.get_name_list():
        try:
            out[name] = cmd.get(name)
        except Exception:
            pass
    return out


def commands():
    """Run core commands and app modules; return their results."""
    results = {}
    cmd.delete('all')
    load('1hpv')
    cmd.ramp_new('ramp', '1hpv', [0, 10], ['blue', 'red'])
    check(results, 'symmetry', lambda: cmd.get_symmetry('1hpv'))
    check(
        results, 'symexp', lambda:
        (cmd.symexp('sym', '1hpv', '1hpv', 5), [(o, cmd.count_atoms(o))
                                                for o in sorted(
                                                    cmd.get_names('objects'))
                                                if o.startswith('sym')])[1])

    def rename_ramp():
        cmd.set_name('ramp', 'ramp2')
        cmd.color('ramp2', '1hpv')
        return cmd.get_names('objects')

    check(results, 'rename ramp', rename_ramp)

    from pymol import xray
    check(results, 'space group P212121',
          lambda: xray.sg_sym_to_mat_list('P212121'))
    check(results, 'space group alias h32', lambda: xray.sg_canonicalize('h32'))

    cmd.delete('all')
    load('1hsg')
    load('1hvr')
    check(results, 'align', lambda: cmd.align('1hvr', '1hsg'))
    check(results, 'distance',
          lambda: cmd.distance('d1', '1hsg///10/CA', '1hsg///20/CA'))
    check(
        results, 'angle',
        lambda: cmd.angle('a1', '1hsg///10/CA', '1hsg///15/CA', '1hsg///20/CA'))
    check(results, 'select', lambda: cmd.select('sele', 'resi 10-20 and 1hsg'))
    check(
        results, 'counts', lambda: {
            s: cmd.count_atoms(s)
            for s in ('all', 'organic', 'solvent', 'polymer', 'hetatm',
                      'name CA', 'sele')
        })
    check(results,
          'color and label',
          lambda:
          (cmd.color('red', 'sele'), cmd.label('sele and name CA', 'resn')),
          keep=False)
    check(
        results, 'scenes', lambda:
        (cmd.scene('s1', 'store'), cmd.set('ambient', 0.25),
         cmd.scene('s2', 'store'), cmd.scene('s1', 'recall', animate=0),
         cmd.get_scene_list())[-1])
    check(
        results, 'movie', lambda: (cmd.mset('1x10'), cmd.mview('store', 1),
                                   cmd.frame(5), cmd.count_frames())[-1])

    def session():
        path = os.path.join(OUT, 'session.pse')
        cmd.save(path)
        cmd.delete('all')
        cmd.load(path)
        return {
            'objects': cmd.get_names('all'),
            'atoms': cmd.count_atoms('all'),
            'scenes': cmd.get_scene_list(),
            'frames': cmd.count_frames(),
        }

    check(results, 'session round trip', session)

    def formats():
        out = {}
        cmd.delete('all')
        load('1hsg')
        cmd.create('ligand', '1hsg and organic')
        for ext in ('pdb', 'sdf', 'mol', 'mol2', 'xyz', 'cif'):
            path = os.path.join(OUT, 'ligand.' + ext)
            cmd.save(path, 'ligand')
            cmd.load(path, 'copy_' + ext)
            out[ext] = [
                cmd.count_atoms('copy_' + ext),
                len(cmd.get_bonds('copy_' + ext)),
                rounded(cmd.get_extent('copy_' + ext), 2),
            ]
        return out

    check(results, 'file format round trips', formats)

    def exports():
        cmd.delete('all')
        load('1hsg')
        cmd.hide('everything')
        cmd.show('cartoon')
        cmd.show('sticks', 'organic')
        cmd.set_view(VIEWS['protease'])
        out = {}
        for name, fn in (('povray', cmd.get_povray), ('vrml1',
                                                      lambda: cmd.get_vrml(1)),
                         ('vrml2', cmd.get_vrml), ('idtf', cmd.get_idtf),
                         ('obj', cmd.get_mtl_obj)):
            text = ''.join(fn())
            out[name] = [len(text), text.count('\n'), zlib.crc32(text.encode())]
        return out

    check(results, 'scene exports', exports)

    def maps():
        cmd.delete('all')
        load('1ubq')
        cmd.map_new('map', 'gaussian', 1.0, '1ubq', 3)
        path = os.path.join(OUT, 'map.ccp4')
        cmd.save(path, 'map')
        cmd.load(path, 'copy')
        out = {'ccp4 round trip': cmd.get_volume_histogram('copy', 16)}
        cmd.map_double('copy')
        out['doubled'] = cmd.get_volume_histogram('copy', 16)
        cmd.map_halve('copy')
        out['halved'] = cmd.get_volume_histogram('copy', 16)
        cmd.map_trim('copy', '1ubq and resi 1-20', 2.0)
        out['trimmed'] = [
            cmd.get_volume_histogram('copy', 16),
            cmd.get_extent('copy')
        ]
        return out

    check(results, 'maps', maps)

    for name in ('mmol', 'mmol.inspector', 'mmol.measure', 'mmol.movie',
                 'mmol.sequence', 'mmol.settings', 'mmol.classic_ray',
                 'mmol.console_hint', 'mmol.move', 'mmol.pick', 'mmol.keys',
                 'mmol.scene_anim', 'mmol.scenes', 'mmol.tmp', 'mmol.tools',
                 'mmol.pymolrc', 'mmol.render_options', 'mmol.mcp',
                 'mmol.mcp.server', 'mmol.mcp.tools', 'mmol.mcp.events',
                 'mmol.mcp.mainthread'):
        check(results,
              'import ' + name,
              lambda name=name: importlib.import_module(name),
              keep=False)

    from mmol import console_hint, inspector, pick, render_options, sequence
    import pymol
    check(results, 'legacy module name in saved movies',
          lambda: pymol.mmol_scene_anim.__name__)
    check(
        results, 'session hooks',
        lambda: sorted(f.__module__
                       for f in pymol._session_save_tasks
                       if f.__module__.startswith('mmol')))
    check(results, 'console hint',
          lambda: console_hint.hint('fetch 1abc, type='))
    check(results, 'mcp tools',
          lambda: len(importlib.import_module('mmol.mcp.tools').TOOLS))
    # These depend on the window size: check only that they run.
    for name, fn in (
        ('pick', lambda: pick.pick_at(0.0, 0.0, 1.5)),
        ('hover preview', lambda: pick.hover_preview_at(0.0, 0.0, 1.5)),
        ('inspector poll',
         lambda: inspector.poll(cmd.get_names(), extent=False)),
        ('sequence poll', sequence.poll),
        ('render options', render_options.emit),
    ):
        check(results, name, fn, keep=False)
    return results


def main():
    """Run the test and write ``report.txt``."""
    report = Report()
    try:
        compare_json(report, 'settings', settings())
        views_path = os.path.join(GOLDEN, 'views.json')
        views = VIEWS
        if not UPDATE:
            views.update(json.load(open(views_path)))
        render(views, 'protease', 400, 300, protease)
        render(views, 'ubiquitin', 500, 360, ubiquitin)
        render(views, 'ramp', 400, 300, ramp)
        render(views, 'dna', 400, 300, dna)
        cmd.space('rgb', quiet=1)
        if UPDATE:
            with open(views_path, 'w') as f:
                json.dump(views, f, indent=1)
                f.write('\n')
        for scene in ('protease', 'ubiquitin', 'ramp', 'dna'):
            for kind in ('gpu', 'cpu'):
                compare_image(report, '%s_%s.png' % (scene, kind))
        compare_json(report, 'commands', commands())
    except Exception:
        report.add(False, 'test script', traceback.format_exc())
    report.write()


main()
