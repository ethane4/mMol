"""Movie export quality overrides (#581).

The native app's Export Movie sheet renders every frame offscreen through the
Metal pipeline. Offline export can afford far more than the 60 fps live view,
so the sheet's Advanced section raises tessellation and anti-aliasing for the
export only. This module owns that override step:

* ``apply_overrides(json)`` snapshots the global value of every setting it is
  about to change, then sets it. The geometry rebuild happens once, on the
  exporter's first on-main ``updateScene``.
* ``restore()`` puts every snapshotted value back. The exporter calls it after
  the export completes, is cancelled, or fails, so the user's session is never
  left changed.

``movie_export`` is the scripting entry point (MCP, the copilot, scripts). It
does not render anything itself: rendering needs the Swift exporter, so it
emits a ``MOVIEEXPORT:<json>`` feedback line that PyMOLEngine picks up and
hands to a MovieExporter. The export runs asynchronously; completion is
reported in the log as `` movie_export: wrote <path>`` (failures as
`` movie_export: failed: <reason>``).

Kept free of Swift-specific state so the preset tables and the snapshot/restore
logic are testable under plain PyMOL.
"""

import json
import os

from pymol import cmd

# Settings the exporter may override. Anything else in a request is rejected,
# so a stray key can't leave an unrelated setting changed after export.
OVERRIDABLE = (
    'cartoon_sampling',
    'cartoon_oval_quality',
    'cartoon_tube_quality',
    'cartoon_loop_quality',
    'ribbon_sampling',
    'sphere_quality',
    'stick_quality',
    'cgo_sphere_quality',
    'surface_quality',
    'metal_msaa',
)

# Quality presets. 'standard' is empty on purpose: no overrides means the
# export is byte-for-byte what it was before the Advanced section existed.
QUALITY_PRESETS = {
    'draft': {
        'metal_msaa': 0,
    },
    'standard': {},
    'high': {
        'cartoon_sampling': 14,
        'cartoon_oval_quality': 16,
        'cartoon_tube_quality': 12,
        'cartoon_loop_quality': 10,
        'ribbon_sampling': 4,
        'sphere_quality': 3,
        'stick_quality': 16,
        'cgo_sphere_quality': 3,
        'metal_msaa': 1,
    },
    'maximum': {
        'cartoon_sampling': 20,
        'cartoon_oval_quality': 24,
        'cartoon_tube_quality': 16,
        'cartoon_loop_quality': 12,
        'ribbon_sampling': 8,
        'sphere_quality': 4,
        'stick_quality': 24,
        'cgo_sphere_quality': 4,
        'surface_quality': 2,
        'metal_msaa': 1,
    },
}

# Supersampling factors per preset (render larger, downsample in the encoder).
PRESET_SUPERSAMPLE = {'draft': 1, 'standard': 1, 'high': 1, 'maximum': 2}

# Global values captured by apply_overrides, keyed by setting name.
_saved = {}


def preset_overrides(quality):
    """The override dict for a named preset (a copy; unknown names raise)."""
    key = str(quality).strip().lower()
    if key not in QUALITY_PRESETS:
        raise ValueError('unknown quality %r (expected one of: %s)' %
                         (quality, ', '.join(QUALITY_PRESETS)))
    return dict(QUALITY_PRESETS[key])


def _validate(overrides):
    bad = [k for k in overrides if k not in OVERRIDABLE]
    if bad:
        raise ValueError('not an export override: %s' % ', '.join(sorted(bad)))


def apply_overrides(overrides):
    """Snapshot then set each override (a dict or a JSON string).

    Safe to call when a previous snapshot is still held: values already saved
    keep their ORIGINAL snapshot, so a second apply can't capture an override
    as the value to restore."""
    if isinstance(overrides, str):
        overrides = json.loads(overrides) if overrides.strip() else {}
    overrides = dict(overrides or {})
    _validate(overrides)
    for name, value in overrides.items():
        if name not in _saved:
            _saved[name] = cmd.get(name)
        cmd.set(name, value, quiet=1)


def restore():
    """Put every snapshotted setting back. Idempotent."""
    errors = []
    for name, value in list(_saved.items()):
        try:
            cmd.set(name, value, quiet=1)
        except Exception as e:  # keep restoring the rest
            errors.append('%s: %s' % (name, e))
        del _saved[name]
    if errors:
        print('MOVIEEXPORT_ERR:restore failed for ' + '; '.join(errors))


def saved_settings():
    """The current snapshot (for tests and diagnostics)."""
    return dict(_saved)


def _emit(line):
    """Send a line to the app's feedback stream. Writes to pcatch (the
    embedded app's stdout->feedback sink) directly when present, so a caller
    that captures stdout -- MCP run_python, the copilot -- can't swallow the
    request. Falls back to print under plain PyMOL."""
    try:
        import pcatch
        pcatch.write(line + '\n')
    except ImportError:
        print(line)


_FORMATS = ('mp4', 'mov', 'gif', 'png')
_VIDEO_CODECS = ('h264', 'hevc', 'prores')


def _infer(filename, format, codec):
    ext = os.path.splitext(filename)[1].lower().lstrip('.')
    fmt = (format or '').lower() or (ext if ext in _FORMATS else 'png')
    if fmt not in _FORMATS:
        raise ValueError('format must be one of: ' + ', '.join(_FORMATS))
    c = (codec or '').lower()
    if fmt == 'gif':
        c = ''
    elif fmt == 'png':
        c = 'png'
    elif not c:
        c = 'prores' if fmt == 'mov' else 'h264'
    if fmt in ('mp4', 'mov') and c not in _VIDEO_CODECS:
        raise ValueError('codec must be one of: ' + ', '.join(_VIDEO_CODECS))
    if c == 'prores' and fmt == 'mp4':
        raise ValueError('ProRes needs a .mov file')
    return fmt, c


def movie_export(filename,
                 width=1920,
                 height=1080,
                 quality='standard',
                 format='',
                 codec='',
                 fps=0,
                 supersample=0,
                 bitrate=0,
                 first=1,
                 last=0,
                 quiet=1,
                 _self=cmd):
    '''"movie_export" renders the movie through mMol's Metal exporter, with export-only quality overrides.

    The live view and session settings are restored when it finishes. Runs asynchronously; the log reports when the file is written. Requires the mMol app (not command-line PyMOL).

    Parameters
    ----------
    filename : str
        output path. The extension picks the format: .mp4
        (H.264 or HEVC), .mov (ProRes 422 on macOS, or H.264/HEVC), .gif, or a
        new or empty folder for a PNG sequence. A movie file replaces an existing
        file; a PNG sequence never writes into a non-empty folder.

        width, height = int: output size in pixels {default: 1920, 1080}
    quality : optional
        draft, standard, high, or maximum {default: standard}
    format : optional
        mp4, mov, gif, or png {default: from the extension}
    codec : optional
        h264, hevc, prores {default: h264 for .mp4, prores for .mov}
    fps : int, optional
        frame rate; 0 uses movie_fps {default: 0}
    supersample : optional
        1, 2 or 4: render larger and downsample; 0 uses the
        quality preset's factor {default: 0}
    bitrate : float, optional
        H.264/HEVC average bitrate in Mbit/s (1-400); 0 lets
        the encoder decide {default: 0}

        first, last = int: frame range; last 0 means the final frame

    Notes
    -----
    Command-line usage::

        movie_export filename [, width [, height [, quality [, format [, codec
            [, fps [, supersample [, bitrate ]]]]]]]]

    Examples
    --------
    movie_export ~/Desktop/spin.mp4, 3840, 2160, quality=high
    movie_export ~/Desktop/spin.mov, 1920, 1080, codec=prores
    '''
    if not str(filename).strip():
        raise ValueError('filename is required')
    filename = os.path.abspath(os.path.expanduser(str(filename)))
    fmt, c = _infer(filename, format, codec)
    q = str(quality).strip().lower()
    overrides = preset_overrides(q)
    ss = int(supersample) or PRESET_SUPERSAMPLE[q]
    if ss not in (1, 2, 4):
        raise ValueError('supersample must be 1, 2 or 4')
    br = float(bitrate)
    if not (br == 0 or 1 <= br <= 400):  # also rejects nan/inf
        raise ValueError('bitrate must be 0 (automatic) or 1-400 Mbit/s')
    first, last = int(first), int(last)
    if first < 1 or last < 0 or (last and last < first):
        raise ValueError('need 1 <= first <= last (last 0 = final frame)')
    w, h = int(width), int(height)
    if w < 16 or h < 16 or w > 8192 or h > 8192:
        raise ValueError('width and height must be between 16 and 8192')
    # The movie's length and rate as the core sees them NOW: the app's own
    # playback mirror lags by a poll, so a script that just ran mset would
    # otherwise look like it has no movie.
    frames = int(_self.count_frames())
    fps = int(fps) or int(round(float(_self.get('movie_fps'))))
    request = {
        'path': filename,
        'frames': frames,
        'width': w,
        'height': h,
        'format': fmt,
        'codec': c,
        'quality': q,
        'overrides': overrides,
        'supersample': ss,
        'fps': fps,
        'bitrate': br,
        'first': first,
        'last': last,
    }
    _emit('MOVIEEXPORT:' + json.dumps(request, separators=(',', ':')))
    if not int(quiet):
        print(' movie_export: queued %dx%d %s (%s) -> %s' %
              (w, h, fmt, q, filename))
    return request
