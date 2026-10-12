'''PyMOL Molecular Graphics System
Copyright (c) Schrodinger, Inc.

The mMol app embeds Python and imports this package; it is not launched
from the command line.
'''

import os
import sys

import _thread as thread

import copy
import threading
import re
import time
import traceback
import math

from . import invocation
from . import colorprinting


def _init_internals(_pymol):

    # Create a temporary object "stored" in the PyMOL global namespace
    # for usage with evaluate based-commands such as alter

    _pymol.stored = Scratch_Storage()

    # Create a permanent object in the PyMOL global namespace
    # that will be picked and unpickled along with the session

    _pymol.session = Session_Storage()

    # This global will be non-None if logging is active
    # (global variable used for efficiency)

    _pymol._log_file = None

    # This global will be non-None if an external gui
    # exists. It mainly exists so that events which occur
    # in the Python thread can be handed off to the
    # external GUI thread through one or more FIFO Queues
    # (global variable used for efficiency)

    # lists of functions to call when saving and restoring pymol session objects

    _pymol._session_save_tasks = []
    _pymol._session_restore_tasks = []

    # cached results (as a list):
    # [ [size, (hash1, hash2, ... ), (inp1, inp2, ...), output],
    #   [size, (hash1, hash2, ... ), (inp1, inp2, ...), output],
    #   ... ]

    _pymol._cache = []

    # stored views

    _pymol._view_dict = {}
    _pymol._view_dict_sc = None

    # stored scenes

    _pymol._scene_quit_on_action = ''

    # get us a private invocation pseudo-module

    _pymol._invocation = Scratch_Storage()
    _pymol._invocation.options = copy.deepcopy(invocation.options)
    _pymol._invocation.get_user_config = invocation.get_user_config

    # these locks are to be shared by all PyMOL instances within a
    # single Python interpeter

    _pymol.lock_api = threading.RLock()  # mutex for API calls from the outside
    _pymol.lock_api_status = threading.RLock()  # mutex for PyMOL status info
    _pymol.lock_api_glut = threading.RLock()  # mutex for GLUT avoidance
    _pymol.lock_api_data = threading.RLock(
    )  # mutex for internal data structures


def get_version_message(v=None):
    '''Get an informative product + version string
    '''
    if not v:
        v = _cmd.get_version()

    p = "PyMOL %s " % v[0]
    p += "Incentive Product" if invocation.options.incentive_product else \
         "Open-Source"

    if v[4]:
        p += ' (' + v[4][:10] + ')'

    if v[3]:
        p += ', ' + time.strftime('%Y-%m-%d', time.localtime(v[3]))

    return p


def exec_str(self, string):
    '''Execute string in "self" namespace (used from C)
    '''
    try:
        exec(string, self.__dict__, self.__dict__)
    except Exception:
        traceback.print_exc()
    return None


def exec_deferred(self):
    '''Run once the core is ready (called from C)
    '''
    self.cmd.config_mouse(quiet=1)


def adapt_to_hardware(self):
    '''Use every CPU core, and store that as the default state (called from C)
    '''
    cmd = self.cmd
    ncpu = os.cpu_count()
    cmd.set("max_threads", ncpu)
    if invocation.options.show_splash:
        print(" Detected %d CPU cores.  Enabled multithreaded rendering." % ncpu)
    cmd.reinitialize("store")


class CmdException(Exception):
    '''Exception type for PyMOL commands
    '''
    label = "Error"

    def __init__(self, message='', label=None):
        self.message = message
        if message:
            self.args = (message,)
        if label:
            self.label = label

    def __str__(self):
        return " %s: %s" % (self.label, self.message)


class IncentiveOnlyException(CmdException):
    '''Exception type for features that are not available in Open-Source PyMOL
    '''
    label = "Incentive-Only-Error"

    def __init__(self, message=''):
        if not message:
            try:
                funcname = sys._getframe(1).f_code.co_name
                message = '"%s" is not available in Open-Source PyMOL' % (
                    funcname,)
            except:
                message = 'Not available in Open-Source PyMOL'
        message += '\n\n' \
                    '    Please visit http://pymol.org if you are interested in the\n' \
                    '    full featured "Incentive PyMOL" version.\n'
        super(IncentiveOnlyException, self).__init__(message)


class Scratch_Storage:
    '''Generic namespace
    '''

    def __reduce__(self):
        # for loading Python 3 (new-style class) pickle with Python 2
        return (self.__class__, (), self.__dict__)

    def get_unused_name(self, prefix='tmp'):
        '''Get an unused name from this namespace
        '''
        i = 1
        while True:
            name = prefix + str(i)
            if not hasattr(self, name):
                setattr(self, name, None)
                return name
            i += 1


class Session_Storage:
    '''Generic namespace
    '''

    def __reduce__(self):
        # for loading Python 3 (new-style class) pickle with Python 2
        return (self.__class__, (), self.__dict__)


def _colortype(cmd):
    # backwards compatible color index type for iterate, which used
    # to expose colors as RGB tuples
    get_color_tuple = cmd.get_color_tuple

    class Color(int):

        def __getitem__(self, i):
            return get_color_tuple(self)[i]

        def __len__(self):
            return 3

    return Color


######### VARIABLES ############################

glutThread = None

# initialize instance-specific module/object internals
_init_internals(sys.modules[__name__])

########## C MODULE ############################

import pymol._cmd

_cmd = sys.modules['pymol._cmd']

get_capabilities = _cmd.get_capabilities

from . import cmd

cmd._COb = None

try:
    import epymol
except ImportError:
    pass

########## WORKAROUND TO PREVENT "import cmd" ##############################
# Previous versions of PyMOL did relative imports and thus allowd
# "import cmd" in pymol scripts to import the pymol.cmd module. To be more
# strict and for compatibility with python3 we use absolute imports now,
# which unfortunately will import an unrelated "cmd" module from the default
# python library, and even worse will corrupt the pymol namespace with it.
# The following causes an import error for "import cmd":


class _NoCmdFinder:

    def find_spec(self, fullname, path=None, target=None):
        if path is None and fullname == 'cmd':
            msg = 'use "from pymol import cmd" instead of "import cmd"'
            print('Warning: {}'.format(msg))
        return None

    find_module = find_spec


sys.meta_path.insert(0, _NoCmdFinder())

########## LEGACY PRINT STATEMENT FOR PYMOL COMMAND LINE ###################

if True:

    def _print_statement(*args, **_):
        '''Legacy Python-2-like print statement for the PyMOL command line'''
        kw = {}
        if args and args[0].startswith('>>'):
            kw['file'] = eval(args[0][2:])
            args = args[1:]
        if args and not args[-1]:
            kw['end'] = ' '
            args = args[:-1]
        args = [eval(a) for a in args]
        print(*args, **kw)

    cmd.extend('print', _print_statement)

# mMol's app modules moved from this package to `mmol`. Movies saved before
# that store frame commands such as "from pymol import mmol_scene_anim", so
# that name still resolves (PEP 562 module __getattr__).
_MOVED_TO_MMOL = {'mmol_scene_anim': 'mmol.scene_anim'}


def __getattr__(name):
    target = _MOVED_TO_MMOL.get(name)
    if target is None:
        raise AttributeError("module 'pymol' has no attribute %r" % name)
    import importlib
    return importlib.import_module(target)
