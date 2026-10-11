"""Console hints for the mMol> prompt, in the style of the fish shell: a flat
signature with the argument being typed, what that argument means, an example,
and completions that carry a short description.

hint(line) -> JSON {
  cmd, about        the command and a one-line summary
  params            [{name, req, default}]  (no `quiet`/`_self`)
  arg               index into params of the argument being typed (-1: none)
  argdoc            what that argument means (and its values)
  example           one example line
  choices           [{text, desc, post}]: completions for the current word,
                    which starts at `start` in line; `post` follows a pick
}

Argument help comes from a short table for the commands people type most
(COMMON), else from the docstring's ARGUMENTS section and the signature.
"""

import glob
import inspect
import json
import os
import re

from pymol import cmd

_lists = re.compile(r"\[[^\]]*\]")
_FILE_CMDS = {
    'load', 'run', 'save', 'cd', 'ls', 'png', 'load_traj', 'spawn', 'mpng'
}
_HIDDEN = {'quiet', '_self', 'kwargs', 'flags'}
MAX = 60

_REPS = 'lines sticks spheres cartoon ribbon surface mesh dots labels nonbonded nb_spheres cell everything'.split(
)
_SEL = 'atoms: an object, chain A, resi 10-20, resn HEM, sele, …'

# Plain-language help for the commands typed most: about, per-argument help
# (text, or (text, [(value, what it is)]) whose values become completions),
# and an example.
_REP_VALUES = [('cartoon', 'secondary-structure cartoon'),
               ('sticks', 'bonds as sticks'), ('lines', 'thin bonds'),
               ('spheres', 'space-filling atoms'),
               ('surface', 'molecular surface'), ('mesh', 'surface as a mesh'),
               ('ribbon', 'backbone trace'), ('dots', 'dotted surface'),
               ('labels', 'atom labels'),
               ('nonbonded', 'unbonded atoms (waters, ions)'),
               ('nb_spheres', 'unbonded atoms as spheres'),
               ('cell', 'unit cell'), ('everything', 'all of the above')]
COMMON = {
    'fetch': ('Download structures (PDB ID, assembly, density map) and load them', {
        'code':
            'PDB ID(s): 1ubq, or "1ubq 4hhb"; 5 characters fetch one chain (1a00A)',
        'name':
            'object name (default: the ID)',
        'state':
            'state to load into (0 = append a new state)',
        'type': ('what to download', [
            ('cif', 'the structure (default)'),
            ('assembly', 'biological assembly: the full complex (= pdb1)'),
            ('density', '2Fo-Fc electron density map (= 2fofc)'),
            ('difference', 'Fo-Fc difference density map (= fofc)'),
            ('em', 'cryo-EM map by EMDB ID, e.g. fetch 1234 (= emd)'),
            ('pubchem', 'PubChem compound by CID (= cid)'),
            ('pdb', 'the structure, legacy PDB format'),
            ('bcif', 'the structure, binary CIF'),
            ('pdb1', 'biological assembly 1 (pdb2, pdb3: others)'),
            ('2fofc', '2Fo-Fc density map'), ('fofc', 'Fo-Fc difference map'),
            ('emd', 'cryo-EM map'), ('cid', 'PubChem compound'),
            ('sid', 'PubChem substance')
        ]),
        'async': ('wait for the download?', [('0', 'wait (mMol default)'),
                                             ('1', 'in the background')]),
        'path':
            'folder to download into'
    }, 'fetch 4hhb, type=assembly'),
    'load': ('Open a file (structure, map, session, …)', {
        'filename': 'the file; ~ works',
        'object': 'object name (default: from the file name)',
        'state': 'state to load into (0 = append)',
        'format': 'file format, if the extension is wrong'
    }, 'load ~/Downloads/model.pdb, mymodel'),
    'color': ('Color atoms', {
        'color': 'a color name (red, skyblue, …) or 0xRRGGBB',
        'selection': _SEL
    }, 'color skyblue, chain A'),
    'show': ('Show a representation', {
        'representation': ('what to draw', _REP_VALUES),
        'selection': _SEL
    }, 'show sticks, resn HEM'),
    'hide': ('Hide a representation', {
        'representation': ('what to hide', _REP_VALUES),
        'selection': _SEL
    }, 'hide everything, solvent'),
    'as': ('Show only this representation', {
        'representation': ('what to draw', _REP_VALUES),
        'selection': _SEL
    }, 'as cartoon'),
    'select': ('Make a named selection', {
        'name': 'name for the selection',
        'selection': 'which atoms: ' + _SEL[7:]
    }, 'select site, byres (resn HEM around 5)'),
    'set': ('Change a setting', {
        'name': 'the setting (Tab lists them)',
        'value': 'its new value',
        'selection': 'object it applies to (default: everything)'
    }, 'set cartoon_transparency, 0.5, 1ubq'),
    'get': ('Print a setting', {
        'name': 'the setting',
        'selection': 'object (default: global)'
    }, 'get ray_trace_mode'),
    'align': ('Superpose by sequence, then structure', {
        'mobile': 'what moves',
        'target': 'what it moves onto',
        'cycles': 'outlier-rejection rounds',
        'object': 'name for an alignment object'
    }, 'align 1ubi, 1ubq'),
    'super': ('Superpose by structure (ignores sequence)', {
        'mobile': 'what moves',
        'target': 'what it moves onto'
    }, 'super 1ubi, 1ubq'),
    'cealign': ('Superpose distant homologs (CE)', {
        'target': 'what stays put',
        'mobile': 'what moves'
    }, 'cealign 1ubq, 1ubi'),
    'ray': ('Ray-trace the view', {
        'width':
            'image width in pixels',
        'height':
            'image height in pixels',
        'antialias': ('smoothing', [('0', 'off'), ('1', 'standard'),
                                    ('2', 'high'), ('3', 'higher'),
                                    ('4', 'highest')])
    }, 'ray 2000, 1500'),
    'png': ('Save the view as a PNG', {
        'filename':
            'where to write it',
        'width':
            'pixels (or 10cm with dpi)',
        'height':
            'pixels',
        'dpi':
            'print resolution',
        'ray': ('ray trace first?', [('1', 'ray-traced image'),
                                     ('0', 'the live view')])
    }, 'png ~/Desktop/figure.png, width=2000, dpi=300, ray=1'),
    'save': ('Write atoms or the session to a file', {
        'filename':
            'file; the extension picks the format (.pdb, .cif, .pse, …)',
        'selection':
            _SEL,
        'state':
            'state to write (-1 = current, 0 = all)'
    }, 'save ~/Desktop/chainA.pdb, chain A'),
    'zoom': ('Fit the view to atoms', {
        'selection': _SEL,
        'buffer': 'extra margin in Å',
        'animate': 'seconds (-1 = default)'
    }, 'zoom resi 20-30, 5'),
    'orient': ('Fit and align the view to atoms', {
        'selection': _SEL,
        'animate': 'seconds (-1 = default)'
    }, 'orient ligand'),
    'center': ('Center the view on atoms', {
        'selection': _SEL
    }, 'center resi 42'),
    'cartoon': ('Cartoon style', {
        'type':
            ('the style', [('automatic', 'helices, sheets and loops (default)'),
                           ('loop', 'thin tube everywhere'),
                           ('rect', 'flat rectangles'), ('oval', 'flat ovals'),
                           ('tube', 'round tube'), ('arrow', 'sheet arrows'),
                           ('dumbbell', 'dumbbell helices'),
                           ('putty', 'thickness by B-factor'),
                           ('skip', 'leave out')]),
        'selection': _SEL
    }, 'cartoon putty'),
    'spectrum': ('Color by a gradient', {
        'expression': ('what varies', [('count', 'position along the chain'),
                                       ('b', 'B-factor'), ('q', 'occupancy'),
                                       ('pc', 'partial charge'),
                                       ('resi', 'residue number')]),
        'palette': ('the colors', [('rainbow', 'blue to red'),
                                   ('blue_white_red', 'diverging'),
                                   ('green_white_magenta', 'diverging'),
                                   ('yellow_red', 'heat'),
                                   ('cyan_white_magenta', 'diverging')]),
        'selection': _SEL
    }, 'spectrum b, blue_white_red'),
    'distance': ('Measure a distance (or polar contacts)', {
        'name':
            'name for the measurement',
        'selection1':
            'first atom(s)',
        'selection2':
            'second atom(s)',
        'cutoff':
            'Å (mode 0)',
        'mode': ('what to measure', [('0', 'any pairs within cutoff'),
                                     ('1', 'closest pairs'),
                                     ('2', 'polar contacts (H-bonds)')])
    }, 'distance d1, resi 10 and name CA, resi 20 and name CA'),
    'label': ('Label atoms', {
        'selection': _SEL,
        'expression': 'text: resn, resi, name, "free text", or a mix'
    }, 'label name CA, resn + resi'),
    'remove': ('Delete atoms', {
        'selection': _SEL
    }, 'remove solvent'),
    'delete': ('Delete objects or selections', {
        'name': 'object/selection name, or all'
    }, 'delete sele'),
    'create': ('Copy atoms into a new object', {
        'name': 'new object',
        'selection': _SEL
    }, 'create ligand, resn HEM'),
    'extract': ('Move atoms into a new object', {
        'name': 'new object',
        'selection': _SEL
    }, 'extract ligand, resn HEM'),
    'bg_color': ('Background color', {
        'color': ('the color', [('white', ''), ('black', ''),
                                ('grey90', 'light grey'), ('grey50', 'grey')])
    }, 'bg_color white'),
    'turn': ('Rotate the camera', {
        'axis': ('axis', [('x', 'horizontal'), ('y', 'vertical'),
                          ('z', 'into the screen')]),
        'angle': 'degrees'
    }, 'turn y, 90'),
    'split_chains': ('One object per chain', {
        'selection': 'object(s) to split',
        'prefix': 'name prefix'
    }, 'split_chains 4hhb'),
}

_doc_cache = {}


def _doc(name, func):
    """(about, {arg: text}, example) from the docstring, COMMON first."""
    if name in _doc_cache:
        return _doc_cache[name]
    about, args, example = '', {}, ''
    doc = inspect.cleandoc(getattr(func, '__doc__', None) or '')
    # numpydoc: a summary line, then sections titled by a line of dashes
    parts = re.split(r'^(\w[\w ]*)\n-{3,}$', doc, flags=re.M)
    about = ' '.join(parts[0].strip().split('\n\n')[0].split()).rstrip('.')
    for head, text in zip(parts[1::2], parts[2::2]):
        if head == 'Parameters':
            for m in re.finditer(r'^(\w+)(?: : [^\n]*)?\n((?:[ \t]+.*\n?)*)',
                                 text.strip('\n') + '\n',
                                 flags=re.M):
                args[m.group(1).rstrip('_')] = ' '.join(m.group(2).split())
        elif head == 'Examples' and not example:
            example = next(
                (l.strip()
                 for l in text.splitlines()
                 if l.strip() and not l.strip().startswith(('#', '::'))), '')
    if name in COMMON:
        a, more, ex = COMMON[name]
        about, example = a, ex
        args = dict(args, **more)
    _doc_cache[name] = (about, args, example)
    return _doc_cache[name]


def _params(func):
    """[{name, req, default}] from the signature (async_ shown as async)."""
    try:
        sig = inspect.signature(func)
    except (TypeError, ValueError):
        return []
    out = []
    for p in sig.parameters.values():
        if p.name in _HIDDEN or p.name.startswith('_') or p.kind > p.KEYWORD_ONLY \
                or p.kind == p.VAR_POSITIONAL:
            continue
        d = p.default
        out.append(
            dict(name=p.name.rstrip('_'),
                 req=d is p.empty,
                 default='' if d is p.empty or d is None else str(d)))
    return out


def _values(text):
    """A "a|b|c" or "x: a, b, c" value list in an argument's help, else []."""
    m = re.search(r'([\w.+-]+(?:\s*\|\s*[\w.+ =-]+)+)', text)
    if m:
        return [
            v.strip().split(' ')[0] for v in m.group(1).split('|') if v.strip()
        ]
    m = re.search(r':\s*([\w-]+(?:\s*,\s*[\w-]+){2,})', text)
    return [v.strip() for v in m.group(1).split(',')] if m else []


def _text(help_):
    return help_[0] if isinstance(help_, tuple) else help_


def _keys(sc):
    sc = sc() if callable(sc) else sc
    return getattr(sc, 'keywords', None) or []


def _about(name):
    func = cmd.keyword.get(name, [None])[0]
    return _doc(name, func)[0] if func else ''


def hint(line):
    seg_start = line.rfind(';') + 1
    seg = line[seg_start:]
    seg_start += len(seg) - len(seg.lstrip())
    seg = seg.lstrip()
    out = dict(cmd='',
               about='',
               params=[],
               arg=-1,
               argdoc='',
               example='',
               choices=[],
               start=len(line))

    word = re.match(r'[^\s,]*', seg).group()
    if len(word) == len(seg):  # still typing the command name
        names = cmd.kwhash.interpret(word, True) if word else None
        names = [names] if isinstance(names, str) else names or []
        names = sorted(n for n in names if not n.startswith('_'))[:MAX]
        out.update(
            start=seg_start,
            choices=[dict(text=n, desc=_about(n), post=' ') for n in names])
        name, rest = word, None
    else:
        name = word if word in cmd.keyword else cmd.kwhash.interpret(word)
        if not isinstance(name, str):
            return json.dumps(out)
        rest = seg[len(word):]
    if name not in cmd.keyword:
        return json.dumps(out)

    func = cmd.keyword[name][0]
    about, argdocs, example = _doc(name, func)
    params = _params(func)
    if name in COMMON:  # rarely used arguments stay out of the signature
        for p in params:
            p['minor'] = not p['req'] and p['name'] not in COMMON[name][1]
    out.update(cmd=name, about=about, params=params, example=example)
    if rest is None:
        return json.dumps(out)

    pieces = _lists.sub('', rest).split(',')
    piece = rest.rsplit(',', 1)[-1]
    count = len(pieces) - 1
    kw = re.match(r'\s*(\w+)\s*=\s*(.*)$', piece)
    cur = kw.group(2) if kw else re.sub(r'.*[\s]', '',
                                        piece)  # the word being typed
    out['start'] = len(line) - len(cur)
    names = [p['name'] for p in params]
    idx = names.index(kw.group(1).rstrip('_')) if kw and kw.group(1).rstrip('_') in names else \
        (count if not kw and count < len(params) else -1)
    out['arg'] = idx
    pname = params[idx]['name'] if idx >= 0 else ''
    text = argdocs.get(pname, '')
    documented = []
    if isinstance(text, tuple):
        text, documented = text
    if idx >= 0 and params[idx]['default'] not in (
            '', None) and 'default' not in text:
        text += ('; ' if text else '') + 'default %s' % params[idx]['default']
    out['argdoc'] = text

    choices = []
    # values: PyMOL's own completion tables for this position, else the
    # argument's documented value list
    table = cmd.auto_arg[idx] if 0 <= idx < len(cmd.auto_arg) else {}
    descs = dict(documented)
    if documented:
        choices = [
            dict(text=v, desc=d, post='')
            for v, d in documented
            if v.startswith(cur)
        ]
    elif name in table:
        sc, _, post = table[name]
        choices = [
            dict(text=k, desc=descs.get(k, ''), post=post if not kw else '')
            for k in _keys(sc)
            if k.startswith(cur) and not k.startswith('_')
        ]
        choices.sort(key=lambda c: c['text'])
    elif name in _FILE_CMDS and idx == 0:
        path = cmd.as_pathstr(cur)
        choices = [
            dict(text=p + ('/' if os.path.isdir(p) else ''), desc='', post='')
            for p in sorted(glob.glob(os.path.expanduser(path) + '*'))
        ]
    else:
        choices = [
            dict(text=v, desc='', post='')
            for v in _values(text)
            if v.startswith(cur)
        ]
    # option names (fish's --options): the remaining optional arguments
    if not kw and (cur or not choices):
        used = {
            re.match(r'\s*(\w+)', p).group(1).rstrip('_')
            for p in pieces[:-1]
            if '=' in p
        }
        for i, p in enumerate(params):
            if i > count or not p['req']:
                if p['name'] not in used and i >= count and p[
                        'name'].startswith(cur) and cur != p['name']:
                    choices.append(
                        dict(text=p['name'] + '=',
                             post='',
                             desc=_text(argdocs.get(p['name'], ''))[:60]))
    minor = {p['name'] + '=' for p in params if p.get('minor')}
    choices.sort(key=lambda c: c['text'] in minor)  # documented options first
    out['choices'] = choices[:MAX]
    return json.dumps(out)
