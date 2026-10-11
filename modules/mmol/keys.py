# mmol/keys.py
#
# mMol delivers cmd.set_key bindings via an NSEvent monitor that consumes a
# key event only when a binding actually fired (mMol#258). That rule gives the
# user's ~/.pymolrc the last word for free — but it also means a binding can
# quietly take over a key the app itself uses as a menu shortcut, and the user
# would just see the menu item stop responding to its key.
#
# This module names that collision out loud, once, right after ~/.pymolrc runs.
#
# Only the CTRL-letter keys can collide: every other mMol menu shortcut
# carries ⌘, and the Swift classifier passes ⌘ events straight through to the
# menus.

# Canonical key token -> the menu command it would shadow.
#
# Mirrors the ⌃-letter entries of AppShortcuts in
# swiftui/PyMOLViewer/Shared/PyMOLApp.swift — a new ⌃ shortcut there must be
# added here too (AppShortcutsTests pins the set on the Swift side).
APP_SHORTCUTS = {
    'CTRL-M': 'Move Objects (Tools menu)',
    'CTRL-E': 'Measure Distances (Tools menu)',
    'CTRL-S': 'Box Select (Tools menu)',
}


def audit_shadowed(_self=None):
    '''Print one warning per user binding that shadows a mMol menu shortcut.

    Returns the list of lines printed, so callers and tests can inspect the
    result. Never raises: this runs on the launch path, and a broken audit must
    not take the app's startup with it.
    '''
    try:
        if _self is None:
            from pymol import cmd as _self

        mappings = getattr(_self, 'key_mappings', None)
        if not isinstance(mappings, dict):
            return []

        lines = []
        for key, label in sorted(APP_SHORTCUTS.items()):
            # An empty mapping is how set_key(key, '') CLEARS a binding.
            if not mappings.get(key):
                continue
            line = (" mMol: %s is bound by your startup script; it now"
                    " overrides the \"%s\" shortcut (the menu item still works"
                    " by click)." % (key, label))
            print(line)
            lines.append(line)
        return lines
    except Exception as e:
        # The docstring's promise, kept: this runs on the launch path, so a bug
        # here must degrade to "no warning" rather than break startup — including
        # when the reporting print is itself what fails (closed or custom stdout).
        try:
            print(' Warning: mMol shortcut audit failed: %s' % (e,))
        except Exception:
            pass
        return []
