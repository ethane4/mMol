"""Display ▸ Render Options: the `ray` settings, edited on a small sample protein.

The Swift panel draws whatever SPEC lists. Settings are PyMOL's own (global),
so every change carries over to the user's scene; the representations and
cartoon style apply to the sample only. `snapshot()` (run once the startup
scripts are in) records the user's defaults; `reset()` returns to them.
"""

import json
import os
import tempfile

from pymol import cmd

SAMPLE = '_render_sample'

# (group, setting, label, kind, extra): slider (min, max), menu [(value, label)],
# toggle, color, vector (a light direction), rep and cartoon (sample only)
SPEC = [
    ('Sample', 'cartoon', 'Cartoon', 'rep', None),
    ('Sample', 'sticks', 'Sticks', 'rep', None),
    ('Sample', 'lines', 'Lines', 'rep', None),
    ('Sample', 'spheres', 'Spheres', 'rep', None),
    ('Sample', 'surface', 'Surface', 'rep', None),
    ('Sample', 'ribbon', 'Ribbon', 'rep', None),
    ('Sample', 'cartoon_style', 'Cartoon style', 'cartoon',
     [(s, s.capitalize()) for s in ('automatic', 'loop', 'oval', 'rectangle',
                                    'tube', 'dumbbell', 'putty')]),
    ('Outline', 'ray_trace_mode', 'Mode', 'menu', [(0, 'Off'),
                                                   (1, 'Color + outline'),
                                                   (2, 'Outline only'),
                                                   (3, 'Posterized + outline')
                                                  ]),
    ('Outline', 'ray_trace_color', 'Outline color', 'color', None),
    ('Outline', 'ray_trace_gain', 'Gain', 'slider', (0.0, 0.5)),
    ('Outline', 'ray_trace_slope_factor', 'Slope factor', 'slider', (0.0, 2.0)),
    ('Outline', 'ray_trace_depth_factor', 'Depth factor', 'slider', (0.0, 0.5)),
    ('Outline', 'ray_trace_disco_factor', 'Discontinuity factor', 'slider',
     (0.0, 2.0)),
    ('Lighting', 'ambient', 'Ambient', 'slider', (0.0, 1.0)),
    ('Lighting', 'direct', 'Direct', 'slider', (0.0, 1.0)),
    ('Lighting', 'reflect', 'Reflect', 'slider', (0.0, 1.0)),
    ('Lighting', 'power', 'Power', 'slider', (0.1, 4.0)),
    ('Lighting', 'reflect_power', 'Reflect power', 'slider', (0.1, 4.0)),
    ('Lighting', 'ray_direct_shade', 'Direct shade', 'slider', (0.0, 1.0)),
    ('Lighting', 'ray_scatter', 'Scatter', 'slider', (0.0, 1.0)),
    ('Lighting', 'ray_legacy_lighting', 'Legacy lighting', 'slider', (0.0,
                                                                      1.0)),
    ('Lighting', 'light_count', 'Lights', 'slider', (1, 10)),
    ('Lighting', 'light', 'Light direction', 'vector', None),
    ('Lighting', 'two_sided_lighting', 'Two-sided lighting', 'menu',
     [(-1, 'Auto'), (0, 'Off'), (1, 'On')]),
    ('Specular', 'specular', 'Specular', 'slider', (0.0, 1.0)),
    ('Specular', 'spec_power', 'Shininess', 'slider', (1.0, 200.0)),
    ('Specular', 'spec_direct', 'Direct specular', 'slider', (0.0, 1.0)),
    ('Specular', 'spec_direct_power', 'Direct shininess', 'slider', (1.0,
                                                                     200.0)),
    ('Specular', 'spec_count', 'Specular lights', 'slider', (-1, 10)),
    ('Shadows', 'ray_shadows', 'Shadows (CPU)', 'toggle', None),
    ('Shadows', 'ray_shadow_decay_factor', 'Decay', 'slider', (0.0, 1.0)),
    ('Shadows', 'ray_shadow_decay_range', 'Decay range', 'slider', (0.0, 10.0)),
    ('Depth', 'depth_cue', 'Depth cue (fog)', 'toggle', None),
    ('Depth', 'fog', 'Fog density', 'slider', (0.0, 1.0)),
    ('Depth', 'fog_start', 'Fog start', 'slider', (0.0, 1.0)),
    ('Background', 'bg_rgb', 'Background', 'color', None),
    ('Background', 'ray_opaque_background', 'Opaque background', 'menu',
     [(-1, 'Auto'), (0, 'Off'), (1, 'On')]),
    ('Quality', 'antialias', 'Antialias', 'menu', [(0, 'Off'),
                                                   (1, '1'), (2, '2'), (3, '3'),
                                                   (4, '4')]),
    ('Transparency', 'transparency', 'Surface', 'slider', (0.0, 1.0)),
    ('Transparency', 'cartoon_transparency', 'Cartoon', 'slider', (0.0, 1.0)),
    ('Transparency', 'stick_transparency', 'Sticks', 'slider', (0.0, 1.0)),
    ('Transparency', 'sphere_transparency', 'Spheres', 'slider', (0.0, 1.0)),
    ('Transparency', 'transparency_mode', 'Mode', 'menu', [(0, '0'), (1, '1'),
                                                           (2, '2'), (3, '3')]),
    ('Transparency', 'ray_transparency_contrast', 'Contrast', 'slider', (1.0,
                                                                         3.0)),
    ('Transparency', 'ray_transparency_specular', 'Specular', 'slider', (0.0,
                                                                         1.0)),
    ('Transparency', 'ray_transparency_oblique', 'Oblique', 'slider', (0.0,
                                                                       1.0)),
    ('Cartoon', 'cartoon_fancy_helices', 'Fancy helices', 'toggle', None),
    ('Cartoon', 'cartoon_cylindrical_helices', 'Cylindrical helices', 'menu',
     [(0, 'Off'), (1, 'On'), (2, 'Smooth')]),
    ('Cartoon', 'cartoon_flat_sheets', 'Flat sheets', 'toggle', None),
    ('Cartoon', 'cartoon_smooth_loops', 'Smooth loops', 'toggle', None),
    ('Cartoon', 'cartoon_loop_radius', 'Loop radius', 'slider', (0.05, 1.0)),
    ('Cartoon', 'cartoon_oval_width', 'Oval width', 'slider', (0.05, 1.0)),
    ('Cartoon', 'cartoon_oval_length', 'Oval length', 'slider', (0.3, 3.0)),
    ('Cartoon', 'cartoon_rect_width', 'Sheet width', 'slider', (0.1, 1.0)),
    ('Cartoon', 'cartoon_rect_length', 'Sheet length', 'slider', (0.5, 3.0)),
]

_SETTINGS = [s for g, s, l, k, x in SPEC if k not in ('rep', 'cartoon')]
_defaults = {}
_saved = None  # (enabled objects, view) while the sample is shown
_sample_view = None
_style = 'automatic'


def snapshot():
    """Record the current settings as the defaults (after ~/.pymolrc)."""
    _defaults.update({s: cmd.get(s) for s in _SETTINGS})


def _value(s):
    t, v = cmd.get_setting_tuple(s)
    v = v[0] if len(v) == 1 else list(v)
    if t == 5:  # colour index
        return list(cmd.get_color_tuple(v)) if v >= 0 else [0.0, 0.0, 0.0]
    if t == 4 and s == 'bg_rgb':
        return list(v)
    return v


def emit():
    """Send the spec and current values to the panel (RENDEROPTS:<json file>)."""
    values = {s: _value(s) for s in _SETTINGS}
    shown = SAMPLE in cmd.get_names('all')
    for g, s, l, k, x in SPEC:
        if k == 'rep':
            values[s] = shown and cmd.count_atoms('%s and rep %s' %
                                                  (SAMPLE, s)) > 0
    values['cartoon_style'] = _style
    spec = [
        dict(group=g, name=s, label=l, kind=k, extra=x)
        for g, s, l, k, x in SPEC
    ]
    # via a file: the console would wrap a line this long
    path = os.path.join(tempfile.gettempdir(), 'mmol_render_options.json')
    with open(path, 'w') as f:
        json.dump({'spec': spec, 'values': values}, f)
    print('RENDEROPTS:' + path)


def set_option(name, value):
    global _style
    kind = next(k for g, s, l, k, x in SPEC if s == name)
    if kind == 'rep':
        (cmd.show if value else cmd.hide)(name, SAMPLE)
    elif kind == 'cartoon':
        _style = value
        cmd.cartoon(value, SAMPLE)
    elif kind == 'color':
        cmd.set(
            name,
            '0x%02x%02x%02x' % tuple(int(round(c * 255)) for c in value[:3]))
    else:
        cmd.set(name, int(value) if isinstance(value, bool) else value)


def reset():
    for s, v in _defaults.items():
        cmd.set(s, v)
    emit()


def begin(path):
    """Show the sample alone (the user's objects are hidden, not changed)."""
    global _saved
    if _saved is None:
        _saved = (cmd.get_names('objects', enabled_only=1), cmd.get_view())
        cmd.disable('all')
    if SAMPLE not in cmd.get_names('all'):
        cmd.load(path, SAMPLE)
        cmd.remove(SAMPLE + ' and solvent')
        cmd.hide('everything', SAMPLE)
        cmd.show('cartoon', SAMPLE)
        cmd.orient(SAMPLE, animate=0)
    elif _sample_view:
        cmd.set_view(_sample_view)
    cmd.enable(SAMPLE)
    emit()


def end(close=False):
    """Back to the user's scene; `close` also drops the sample."""
    global _saved, _sample_view
    if _saved is not None:
        _sample_view = cmd.get_view()
        cmd.disable(SAMPLE)
        for o in _saved[0]:
            cmd.enable(o)
        cmd.set_view(_saved[1])
        _saved = None
    if close:
        cmd.delete(SAMPLE)
        _sample_view = None
