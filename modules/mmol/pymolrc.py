# mmol/pymolrc.py
#
# mMol boots Python directly and never goes through
# pymol.invocation.parse_args, so ~/.pymolrc would otherwise be silently
# ignored (see mMol#225). load() runs it the way desktop PyMOL does, so one
# startup script serves both apps.

import os

_HOME = os.path.expanduser('~')

PYMOLRC_PY = os.path.join(_HOME, '.pymolrc.py')
PYMOLRC_PML = os.path.join(_HOME, '.pymolrc')


def _seed_main_namespace(_self):
    # The vanilla `pymol` CLI executable IS the process's real __main__
    # module, so it already has `cmd`/`stored` etc. in scope when it runs a
    # user's .pymolrc.py via `run(path, 'main')` — that's what lets pymolrc
    # scripts call bare cmd.xxx() without importing it. mMol boots Python
    # directly, so __main__ starts empty; without this, every .pymolrc.py using
    # that (extremely common) convention would fail with
    # NameError('cmd' is not defined).
    import __main__
    import pymol
    __main__.__dict__.setdefault('cmd', _self)
    __main__.__dict__.setdefault('stored', pymol.stored)


def load(_self=None):
    '''Run ~/.pymolrc.py and ~/.pymolrc, if present.'''
    if _self is None:
        from pymol import cmd as _self

    if os.path.exists(PYMOLRC_PY):
        _seed_main_namespace(_self)
        _self.run(PYMOLRC_PY, 'main')

    if os.path.exists(PYMOLRC_PML):
        _self.do('@' + PYMOLRC_PML)
