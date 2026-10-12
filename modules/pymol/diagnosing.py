import sys
import os


def diagnostics(filename='', compact=0, quiet=1):
    '''Get system level diagnostics data

    Parameters
    ----------
    filename : str, optional
        If given, write output to text file

    Notes
    -----
    Command-line usage::

        diagnostics [ filename ]
    '''
    import time
    from pymol import invocation
    from pymol import cmd, CmdException

    compact, quiet = int(compact), int(quiet)

    TZ = '%+05d' % (time.timezone / 36)

    version = cmd.get_version()
    body = u'PyMOL %s\n' % version[0]

    if not compact:
        if version[3]:
            body += 'build date: %s %s\n' % (time.ctime(version[3]), TZ)
        if version[4]:
            body += 'git sha: %s\n' % version[4]

    body += diagnostics_conda(version[0])

    if not compact:
        body += '\nLicense Information:\n'

    body += 'Open-Source Build\n'

    body += diagnostics_platform()

    if not compact:
        body += '\nOpenGL Driver:\n'
    renderer = cmd.get_renderer()
    body += str(renderer[0] or '(none)') + '\n'
    body += str(renderer[1] or '(none)') + '\n'
    body += str(renderer[2] or '(none)') + '\n'

    if not compact:
        body += diagnostics_python()

        body += '\nStartup Scripts:\n'
        pymolrc = invocation.get_user_config()
        if pymolrc:
            body += '\n'.join(pymolrc) + '\n'
        else:
            body += '(no pymolrc file found)\n'

        body += diagnostics_env_vars()

        body += '\nDiagnostics collected on %s %s\n' % (time.ctime(), TZ)

    if filename:
        filename = cmd.exp_path(filename)
        if not filename.endswith('.txt'):
            raise CmdException('filename must have .txt extension')

        with open(filename, 'w') as handle:
            handle.write(body)
        print('Diagnostics written to "%s"' % filename)
    elif not quiet:
        print(body.rstrip())

    return body


def diagnostics_conda(versionselector=''):
    import glob
    body = ''
    try:
        # pymol conda package version
        condameta = glob.glob(
            os.path.join(sys.prefix, 'conda-meta',
                         'pymol-' + versionselector + '*.json'))
        for jsonfile in condameta:
            import json
            d = json.load(open(jsonfile))
            if not versionselector:
                body += os.path.basename(jsonfile) + '\n'
            body += u'conda build: {} {}\n'.format(
                d.get('build'),
                d.get('schannel') or d.get('channel'))
    except BaseException as e:
        print(e)
    return body


def diagnostics_platform(compact=False):
    import platform
    body = u''

    if not compact:
        body += '\nOperating System:\n'

    body += platform.platform() + '\n'

    return body


def diagnostics_python():
    tmpl = (u'\nPython:\n'
            u'{v}\n'
            u'prefix={p}\n'
            u'executable={exe}\n'
            u'filesystemencoding={enc}\n')
    return tmpl.format(v=sys.version,
                       p=sys.prefix,
                       exe=sys.executable,
                       enc=sys.getfilesystemencoding())


def diagnostics_env_vars():
    import re
    import textwrap

    body = '\nQt, Python and PyMOL Environment Variables:\n'
    for key in sorted(os.environ):
        if re.search(r'^PY|QT|^LANG', key) and not (key == 'PYMOL_SCRIPTS'):
            body += u'%s=%s\n' % (key, os.environ[key])

    body += '\nPATH:\n'
    body += textwrap.fill(os.getenv('PATH', ''), 78) + '\n'
    return body
