#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------

# invocation.py
#
# Startup options. The app sets its own options before PInit; the C core
# reads these attributes when a pymol2 instance is created (PConvertOptions).
#

import re
import os
import glob

pymolrc_pat1 = '.pymolrc*'
pymolrc_pat2 = 'pymolrc*'

ros_pat = 'run_on_startup*'


class generic:
    pass


options = generic()

options.no_gui = 0
options.internal_gui = 1
options.internal_feedback = 1
options.external_gui = 1
options.force_stereo = -1
options.game_mode = 0
options.show_splash = 1
options.read_stdin = 0
options.win_x = 640
options.win_y = 480
options.win_px = 4
options.win_py = 214
options.sigint_handler = 1  # terminate on Ctrl-C?
options.reuse_helper = 0
options.auto_reinitialize = 0
options.keep_thread_alive = 0
options.after_load_script = ""
options.quiet = 0
options.multisample = 0
options.incentive_product = 0
options.window_visible = 1
options.presentation = 0
options.defer_builds_mode = 0
options.full_screen = 0
options.sphere_mode = -1
options.stereo_capable = 0
options.stereo_mode = 0
options.zoom_mode = -1
options.no_quit = 0
options.exit_on_error = 0
options.launch_status = 0
options.gldebug = 0
options.openvr_stub = False
options.blue_line = 0
options.security = 1

script_re = re.compile(r"pymolrc$|\.pml$|\.PML$|\.p1m$|\.P1M$")
py_re = re.compile(r"\.py$|\.pym$|\.PY$|\.PYM$")


def get_user_config():
    for d in [
            os.getcwd(), '$HOME', '$PYMOL_PATH',
            '{}/pymol'.format(
                os.environ.get('XDG_CONFIG_HOME', '$HOME/.config'))
    ]:
        d = os.path.expandvars(d)
        for pat in [pymolrc_pat1, pymolrc_pat2]:
            lst = glob.glob(d + os.sep + pat)
            if lst:
                break
        if lst:
            break
    # global run_on_startup script (not overridden by pymolrc files)
    if "PYMOL_PATH" in os.environ:
        first = glob.glob(os.environ['PYMOL_PATH'] + "/" + ros_pat)
    else:
        first = []
    second = []
    for a in lst:
        if py_re.search(a):
            first.append(a)  # preceeding "_ " cloaks
        elif script_re.search(a):
            second.append(a)  # preceeding "_ " cloaks

    first.sort()
    second.sort()
    return first + second
