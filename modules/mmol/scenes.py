"""Per-scene "extras" snapshot for mMol (settings + object TTT).

Classic PyMOL scenes store the camera + representations + colors but NOT setting
values or per-object Move-mode transforms, so the render "look" and object
arrangement were not captured by `scene ... store`. This module snapshots both
per scene name and re-applies them on recall, persisting them in the .pse via
registered session save/restore tasks (see cmd._deferred_init_pymol_internals):

  * render "look" settings (CAPTURE below) — metal_* / lighting / fog
  * per-object TTT matrices (Move mode) — _scene_ttt

Camera lens / zoom / orthographic / FOV / clip slab are already restored by the
scene's saved view, so they're intentionally NOT captured here (the view owns them).

A partial session restore (`load x.pse, partial=1`) leaves every extra alone:
the core keeps the live scenes then, so the extras keep matching them.

Driven from `cmd.scene` via the central on_scene_action hook (snapshot on
store/update, apply on recall/prev/next, prune on delete, clear_all on clear,
rename on rename); no per-UI-call-site pairing is needed.
"""
from pymol import cmd

# Render "look" settings a scene captures. All are get/set-able globals.
CAPTURE = [
    "metal_ssao",
    "metal_outline",
    "metal_outline_width",
    "metal_msaa",
    "depth_cue",
    "fog",
    "surface_quality",
    "ambient",
    "direct",
    "reflect",
    "specular",
    "shininess",
    "ray_opaque_background",
]

# {scene_name: {setting: value}} — persisted into the .pse via session tasks.
_scene_settings = {}

# Identity TTT (16-float) used when an object has no transform yet / to reset one.
_IDENTITY_TTT = [
    1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0
]

# {scene_name: {obj_name: [16 floats] | None}} — per-object Move-mode TTT.
# None means the object had no TTT (unmoved) when the scene was stored; on recall
# that resets the object to identity. Every object present at store time is a key,
# so recall can reset as well as move. Persisted into the .pse alongside settings.
_scene_ttt = {}

# Reentrancy guard: internal temp-scene machinery (.pse legacy convert, multi-scene
# export) increments this so the cmd.scene hook does NOT capture/apply during its
# throwaway store/recall/clear (see viewing.session_restore_scenes / exporting).
_suspend = 0


class _Suspend:

    def __enter__(self):
        global _suspend
        _suspend += 1
        return self

    def __exit__(self, *exc):
        global _suspend
        _suspend = max(0, _suspend - 1)
        return False


def suspended():
    """Context manager: pause the cmd.scene capture/apply hook for internal
    temp-scene operations."""
    return _Suspend()


# A key no `cmd.scene` name can be (scene names cannot hold NUL), used only
# while preserved() is active.
_LIVE = '\x00live'


class _Preserved:
    """The snapshot lives on the instance and is installed under _LIVE only
    while it is being put back, so nested preserved() blocks each restore
    their own state, and nothing that walks the scene dicts (session_save,
    prune) sees a phantom entry while the block body runs -- only during the
    brief restore in __exit__, which session_save filters anyway."""

    def __init__(self, _self):
        self._self = _self

    def __enter__(self):
        self._settings = _capture(self._self)
        return self

    def __exit__(self, *exc):
        _scene_settings[_LIVE] = self._settings
        try:
            apply_settings(_LIVE, self._self)
        finally:
            _scene_settings.pop(_LIVE, None)
            for key in [k for k in _reported if k[0] == _LIVE]:
                _reported.discard(key)
        return False


def preserved(_self=cmd):
    """Context manager: put every captured setting back as it was when the
    block exits (#508). Nests.

    For code that has to DISPLAY movie frames only to read something off them
    -- appkit_movie._scene_keyframes scrubs every frame to find the scene cuts
    -- because displaying a frame runs its authored commands: each scene
    keyframe's enter_scene and each interpolated `set`. Without this, building
    a movie left the live session carrying whatever the last scrubbed frame
    applied. Restored through the same conditional writes as a recall, so an
    unchanged value costs nothing."""
    return _Preserved(_self)


def _current(_self=cmd):
    try:
        return _self.get("scene_current_name") or ""
    except Exception:
        return ""


def _capture(_self=cmd):
    out = {}
    for s in CAPTURE:
        try:
            out[s] = _self.get(s)
        except Exception:
            pass  # setting absent in this build — skip
    return out


def _objects(_self=cmd):
    """Current objects, WITHOUT groups.

    `get_names('objects')` includes group objects, and `set_object_ttt` expands
    a group name to its members, so capturing a group and then resetting it to
    identity on recall wiped every member's stored Move-mode transform."""
    try:
        objs = _self.get_names('objects') or []
    except Exception:
        return []
    out = []
    for o in objs:
        try:
            if _self.get_type(o) == 'object:group':
                continue
        except Exception:
            pass
        out.append(o)
    return out


# (scene, setting) pairs already reported by apply_settings. A failure that
# repeats every frame is reported once, not once per frame -- and the entries
# are dropped with the scenes they name, because scene names like "A" or "001"
# collide across sessions constantly: a stale entry would silence a REAL first
# failure in a different .pse for the life of the process.
_reported = set()


def _forget_reports(name):
    for key in [k for k in _reported if k[0] == name]:
        _reported.discard(key)


def apply_settings(name, _self=cmd):
    """Make scene `name`'s captured GLOBAL settings current, skipping every
    write that would not change anything.

    Public because the movie animator replays the same payload at each scene
    keyframe (mmol_scene_anim.enter_scene). It used to have its own copy of
    this loop; the two then have to be fixed twice, and the second copy was
    missed exactly once -- see below for why that is expensive.

    Only what CHANGES is written. surface_quality is cRepInvRep -- a full
    surface RE-TESSELLATION -- so an unconditional re-write costs a rebuild on
    every recall even when nothing differs: 127 s per no-op recall on 1aon
    (58 870 atoms) with a surface shown."""
    d = _scene_settings.get(name)
    if not d:
        return
    known = set(CAPTURE)
    for s, v in d.items():
        # A name this build does not capture comes from a .pse written by a
        # different one. Skipped in silence and left in the payload untouched, so
        # re-saving here does not strip it -- warning about it would put a line
        # on the console at every scene keyframe of every pass of a movie.
        if s not in known:
            continue
        try:
            if _setting_differs(s, v, _self):
                _self.set(s, v)
        except Exception as exc:
            # A name this build DOES capture but cannot accept the value for.
            # Worth saying, but once: apply_settings runs at every scene
            # keyframe, so an unconditional print here is a console line per
            # setting per frame for the life of the movie.
            key = (name, s)
            if key in _reported:
                continue
            _reported.add(key)
            try:
                from pymol import colorprinting
                colorprinting.warning(
                    ' scene: %s could not be restored for scene "%s": %s' %
                    (s, name, exc))
            except Exception:
                pass


def _setting_differs(name, value, _self=cmd):
    """True when writing `value` would actually change `name`.

    Compared as the text `cmd.get` returns, the same accessor the capture used.
    Unknown or unreadable settings answer True, which degrades to the old
    unconditional write rather than silently skipping one."""
    try:
        return _self.get(name) != value
    except Exception:
        return True


def _capture_ttt(_self=cmd):
    """Per-object TTT for every current object (None where unmoved); groups
    are excluded (see _objects)."""
    out = {}
    objs = _objects(_self)
    for o in objs:
        try:
            out[o] = _self.get_object_ttt(o)  # list[16] or None
        except Exception:
            out[o] = None
    return out


def _apply_ttt(name, _self=cmd):
    """Restore per-object TTT for scene `name`; reset unmoved objects to identity;
    skip objects that no longer exist."""
    d = _scene_ttt.get(name)
    if not d:
        return
    # Groups filtered here too, not only at capture: a .pse written before this
    # guard existed still has them in its map.
    live = set(_objects(_self))
    for o, ttt in d.items():
        if o not in live:
            continue
        try:
            _self.set_object_ttt(o, list(ttt) if ttt else list(_IDENTITY_TTT))
        except Exception:
            pass


def scene_settings_map(name):
    """Copy of the render settings captured for scene `name` ({} if none)."""
    return dict(_scene_settings.get(name, {}))


def emit_object_motion(name, frame, _self=cmd):
    """Author per-object Move-mode TTT keyframes for scene `name` at movie `frame`:
    apply each captured object TTT, then store an object-matrix mview keyframe
    (freeze=1 — the caller interpolates per object once at the end). Returns the
    list of objects keyframed. [] if no captured TTT / objects absent."""
    d = _scene_ttt.get(name)
    if not d:
        return []
    try:
        live = set(_self.get_names('objects') or [])
    except Exception:
        live = set()
    done = []
    for o, ttt in d.items():
        if o not in live:
            continue
        try:
            _self.set_object_ttt(o, list(ttt) if ttt else list(_IDENTITY_TTT))
            _self.mview('store', object=o, first=int(frame), freeze=1)
            done.append(o)
        except Exception:
            pass
    return done


def restoring_partial(_self=cmd):
    """True while importing.set_session runs the restore tasks of a partial
    restore (`load x.pse, partial=1`). Read defensively: a stand-in cmd has
    no _pymol."""
    try:
        return bool(
            getattr(getattr(_self, '_pymol', None), '_session_partial_restore',
                    False))
    except Exception:
        return False


def snapshot_current(_self=cmd):
    """Capture the current render settings AND per-object TTT for the current
    scene. Call right after `scene ..., store` / `update`."""
    name = _current(_self)
    if name:
        _scene_settings[name] = _capture(_self)
        _scene_ttt[name] = _capture_ttt(_self)
    return name


def apply(name, _self=cmd):
    """Re-apply scene `name`'s captured render settings and per-object TTT."""
    apply_settings(name, _self)
    _apply_ttt(name, _self)


def apply_current(_self=cmd):
    """Apply the now-current scene's settings (after prev/next navigation)."""
    apply(_current(_self), _self)


def prune(_self=cmd):
    """Drop snapshots for scenes that no longer exist (call after a delete)."""
    try:
        live = set(_self.get_scene_list() or [])
    except Exception:
        return
    for name in list(_scene_settings.keys()):
        if name not in live:
            _scene_settings.pop(name, None)
    for name in {k[0] for k in _reported}:
        if name not in live:
            _forget_reports(name)
    for name in list(_scene_ttt.keys()):
        if name not in live:
            _scene_ttt.pop(name, None)


def clear_all(_self=cmd):
    """Forget all snapshots (call after `scene *, clear`)."""
    _scene_settings.clear()
    _scene_ttt.clear()
    _reported.clear()


def rename(old, new, _self=cmd):
    """Re-key snapshots when a scene is renamed (old -> new)."""
    if not old or not new or old == new:
        return
    if old in _scene_settings:
        _scene_settings[new] = _scene_settings.pop(old)
    for key in [k for k in _reported if k[0] == old]:
        _reported.discard(key)
        _reported.add((new, key[1]))
    if old in _scene_ttt:
        _scene_ttt[new] = _scene_ttt.pop(old)


def on_scene_action(key, action, new_key=None, _self=cmd):
    """Central hook called from cmd.scene AFTER the native op completes. `action`
    is already normalized by cmd.scene (update/append -> store, clear -> delete,
    auto-recall -> next). No-op while suspended()."""
    if _suspend:
        return
    if action in ('store', 'insert_after', 'insert_before'):
        snapshot_current(_self)
    elif action in ('recall', 'next', 'previous'):
        apply_current(_self)
    elif action == 'delete':
        if key == '*':
            clear_all(_self)
        else:
            prune(_self)
    elif action == 'rename':
        rename(key, new_key, _self)


# --- .pse persistence (registered in cmd._deferred_init_pymol_internals) ---
def session_save(session, *, _self=cmd):
    session["mmol_scene_settings"] = dict(_scene_settings)
    session["mmol_scene_ttt"] = {k: dict(v) for k, v in _scene_ttt.items()}
    # preserved() only installs its private entry while restoring; a save can
    # never see it, but a session must never carry it.
    session["mmol_scene_settings"].pop(_LIVE, None)
    return 1


def session_restore(session, *, _self=cmd):
    # A partial restore keeps the live scenes (the core skips the session's),
    # so their extras must stay too: clearing here wiped the settings and TTT
    # of every scene that survives the load.
    if restoring_partial(_self):
        return 1
    _scene_settings.clear()
    _scene_ttt.clear()
    # A different .pse can hold a scene with the same NAME and a different bad
    # value; its first failure has to be reported.
    _reported.clear()
    # sessions saved before the mMol rename used raymol_* keys
    d = session.get("mmol_scene_settings", session.get("raymol_scene_settings"))
    if isinstance(d, dict):
        _scene_settings.update(d)
    t = session.get("mmol_scene_ttt", session.get("raymol_scene_ttt"))
    if isinstance(t, dict):
        _scene_ttt.update({k: dict(v) for k, v in t.items()})
    return 1
