"""General undo and redo for mMol, by session snapshots.

PyMOL's own ``undo`` only reverts an edited object's coordinates. This module
records the whole scene instead. Commands that change it are wrapped once at
startup (``install``); a change only marks the scene changed, and the app's
next background poll (``checkpoint``) saves a session snapshot, after the
screen has updated, once changes have stopped for a moment, so a script, a
slider drag or a gizmo drag becomes one step. ``undo`` and ``redo`` restore
the neighbouring snapshot but keep the current camera, so moving the view is
never undone. Selections and mMol's internal objects (names starting with
``_``) are not recorded.
"""

import contextlib
import functools
import pickle
import threading
import time

from pymol import cmd

MAX_STEPS = 20
MAX_BYTES = 200 * 1024 * 1024
# Changes closer together than this become one step.
SETTLE_SECONDS = 0.4

# Commands that change the scene, not just the camera or a selection.
CHANGING = '''
    color set_color spectrum cartoon show hide show_as toggle set unset set_bond
    unset_bond label remove delete create extract copy load load_traj fetch
    alter alter_state h_add h_fill remove_picked bond unbond fuse attach replace
    sort set_name group ungroup order symexp align super cealign fit pair_fit
    intra_fit translate rotate transform_selection transform_object matrix_copy
    matrix_reset enable disable mset madd minsert mdelete mview scene map_new
    isomesh isosurface isodot volume ramp_new gradient slice_new sculpt_iterate
    clean dss protect deprotect mask unmask flag smooth'''.split()

# Commands that set settings for the mouse, not the scene: never recorded.
NEUTRAL = 'config_mouse mouse edit_mode button'.split()

# Commands that only move the camera unless told to move atoms (camera=0).
_CAMERA_BY_DEFAULT = {'translate': 3, 'rotate': 4}

_lock = threading.RLock()
_local = threading.local()
_history = []  # [(label, pickled session)], oldest first
_pos = -1  # index of the snapshot that matches the scene
_pending = []  # labels of changes not yet saved
_last_change = 0.0
_quiet = False  # while saving or restoring a snapshot


def _snapshot():
    """Return the current session, pickled."""
    return pickle.dumps(cmd.get_session(), protocol=pickle.HIGHEST_PROTOCOL)


def _push(label, blob):
    """Make `blob` the newest step, dropping redo steps and old steps."""
    global _pos
    del _history[_pos + 1:]
    _history.append((label, blob))
    _pos = len(_history) - 1
    while len(_history) > MAX_STEPS + 1 or (len(_history) > 2 and sum(
            len(b) for _, b in _history) > MAX_BYTES):
        _history.pop(0)
        _pos -= 1


def _save_pending():
    """Save the changes made since the last snapshot as one step."""
    global _pending, _quiet
    if not _pending:
        return
    label = _pending[0] + (' …' if len(_pending) > 1 else '')
    _pending = []
    _quiet = True
    try:
        _push(label, _snapshot())
    finally:
        _quiet = False


def _records(name, args, kwargs):
    """Whether a call to command `name` is an undoable change.

    Parameters
    ----------
    name : str
        The command's name.
    args : tuple
        Its positional arguments.
    kwargs : dict
        Its keyword arguments.

    Returns
    -------
    bool
        False for camera-only moves, selections, and mMol's internal objects.
    """
    if name in _CAMERA_BY_DEFAULT:
        i = _CAMERA_BY_DEFAULT[name]
        camera = kwargs.get('camera', args[i] if len(args) > i else 1)
        if int(camera):
            return False
    strings = [
        a for a in list(args) + list(kwargs.values()) if isinstance(a, str)
    ]
    if any(s.startswith('_') for s in strings):
        return False
    if name in ('enable', 'disable', 'delete', 'set_name') and strings:
        try:
            if cmd.get_type(strings[0]) == 'selection':
                return False
        except Exception:
            pass
    return True


def _label(name, args, kwargs):
    """The command as the console would show it, shortened."""
    parts = [str(a) for a in args[:3]]
    parts += ['%s=%s' % (k, v) for k, v in kwargs.items() if k != '_self']
    text = '%s %s' % (name, ', '.join(parts)) if parts else name
    return text if len(text) <= 60 else text[:59] + '…'


def _wrap(fn, name):
    """Return `fn` recording its top-level calls as undoable changes."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        global _quiet, _last_change
        depth = getattr(_local, 'depth', 0)
        top = depth == 0 and not _quiet and _records(name, args, kwargs)
        if top:
            with _lock:
                if not _history:  # the scene before the first change
                    _quiet = True
                    try:
                        _push('start', _snapshot())
                    finally:
                        _quiet = False
        _local.depth = depth + 1
        try:
            return fn(*args, **kwargs)
        finally:
            _local.depth = depth
            if top:
                with _lock:
                    _pending.append(_label(name, args, kwargs))
                    _last_change = time.monotonic()

    wrapper._mmol_undo = True
    return wrapper


def _unrecorded(fn):
    """Return `fn` running with recording suspended."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with suspended():
            return fn(*args, **kwargs)

    wrapper._mmol_undo = True
    return wrapper


def checkpoint():
    """Save pending changes once they have settled (called by the UI poll)."""
    with _lock:
        if _pending and time.monotonic() - _last_change >= SETTLE_SECONDS:
            _save_pending()


@contextlib.contextmanager
def suspended():
    """Run a block whose changes are not recorded (temporary UI state)."""
    global _quiet
    with _lock:
        was, _quiet = _quiet, True
    try:
        yield
    finally:
        _quiet = was


def _restore(index):
    """Restore snapshot `index`, keeping the current camera."""
    global _quiet
    view = cmd.get_view()
    _quiet = True
    try:
        cmd.set_session(pickle.loads(_history[index][1]))
        cmd.set_view(view)
    finally:
        _quiet = False


def undo(_self=cmd):
    """Undo the last change to the scene (the ``undo`` command)."""
    global _pos
    with _lock:
        _save_pending()
        if _pos <= 0:
            print(' Undo: nothing to undo.')
            return
        label = _history[_pos][0]
        _restore(_pos - 1)
        _pos -= 1
    print(' Undo: %s' % label)


def redo(_self=cmd):
    """Redo the last undone change (the ``redo`` command)."""
    global _pos
    with _lock:
        _save_pending()
        if _pos >= len(_history) - 1:
            print(' Redo: nothing to redo.')
            return
        _restore(_pos + 1)
        _pos += 1
        label = _history[_pos][0]
    print(' Redo: %s' % label)


def install():
    """Wrap the changing commands and make ``undo``/``redo`` use this module.

    Safe to call more than once.
    """
    replaced = {}
    for name in CHANGING:
        fn = getattr(cmd, name, None)
        if fn is None or getattr(fn, '_mmol_undo', False):
            continue
        replaced[fn] = _wrap(fn, name)
        setattr(cmd, name, replaced[fn])
    for name in NEUTRAL:
        fn = getattr(cmd, name, None)
        if fn is None or getattr(fn, '_mmol_undo', False):
            continue
        replaced[fn] = _unrecorded(fn)
        setattr(cmd, name, replaced[fn])
    for name, fn in (('undo', undo), ('redo', redo)):
        replaced[getattr(cmd, name)] = fn
        setattr(cmd, name, fn)
    # The command parser looks commands up in cmd.keyword, which holds the
    # functions themselves, including aliases (``as`` for show_as).
    for entry in cmd.keyword.values():
        if entry and entry[0] in replaced:
            entry[0] = replaced[entry[0]]
    # PyMOL tunes settings to the hardware on its first idle pass, which can
    # come after this; that is setup, not a change to undo.
    import pymol
    if not getattr(pymol.adapt_to_hardware, '_mmol_undo', False):
        pymol.adapt_to_hardware = _unrecorded(pymol.adapt_to_hardware)
