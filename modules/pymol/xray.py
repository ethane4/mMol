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

# xray.py
# This section contains python code for supporting
# x-ray crystallography functions


def sg_canonicalize(sg):
    import re
    sg = re.sub(r'\s+', ' ', sg.strip().upper())
    return space_group_map.get(sg, sg)


def sg_register_if_unknown(sg, sym_op):
    sg = sg_canonicalize(sg)
    if sg not in sym_dict:
        sym_op = [op.replace(' ', '') for op in sym_op]
        sym_dict[sg] = sym_op


def sg_sym_to_mat_list(sgsymbol):  # TODO _self
    result = None
    sym_op = sym_dict.get(sg_canonicalize(sgsymbol))
    if sym_op is not None:
        result = []
        for op in sym_op:
            mat = []
            for expr in op.split(','):
                mat.append(expr_to_vect(expr))
            mat.append([0.0, 0.0, 0.0, 1.0])
            result.append(mat)
    else:
        try:
            from .cmd import _feedback, fb_module, fb_mask
            if (_feedback(fb_module.symmetry, fb_mask.errors)):
                print("Symmetry-Error: Unrecognized space group symbol '" +
                      sgsymbol + "'.")
        except:
            pass
    return result


hex_to_rhom_xHM = {
    #  extended Hermann Mauguin symbol translation
    #  from CCP4's syminfo.lib
    'H 3': 'R 3 :H',
    'H -3': 'R -3 :H',
    'H 3 2': 'R 3 2 :H',
    'H 3 m': 'R 3 m :H',
    'H 3 c': 'R 3 c :H',
    'H -3 2/m': 'R -3 m :H',
    'H -3 m': 'R -3 m :H',
    'H -3 2/c': 'R -3 c :H',
    'H -3 c': 'R -3 c :H',
    'H 3 M': 'R 3 m :H',
    'H 3 C': 'R 3 c :H',
    'H -3 2/M': 'R -3 m :H',
    'H -3 M': 'R -3 m :H',
    'H -3 2/C': 'R -3 c :H',
    'H -3 C': 'R -3 c :H',
}


def expr_to_vect(e):
    '''Example:
    >>> expr_to_vect('x-z+1/2')
    [1., 0., -1., 0.5]
    '''
    import re
    r = [0., 0., 0., 0.]
    for sign, var in re.findall(r'(-?)([^-+]+)', e):
        idx = {'x': 0, 'y': 1, 'z': 2}.get(var, 3)
        if idx == 3:
            r[3] = eval(var, {}, {})
        else:
            r[idx] = 1.
        if sign == '-':
            r[idx] *= -1.
    return r


def _load_tables():
    """Load the space-group tables from ``$PYMOL_DATA/pymol/space_groups.json``.

    Returns
    -------
    dict
        ``sym_base``: ``[ops, names]`` pairs (the symmetry operators of a
        space group and its symbols); ``space_group_map``: short-name aliases.
    """
    import json
    import os
    data_dir = os.environ.get('PYMOL_DATA') or os.path.join(
        os.path.dirname(__file__), '..', '..', 'data')
    with open(os.path.join(data_dir, 'pymol', 'space_groups.json')) as f:
        return json.load(f)


_tables = _load_tables()

# Symbol -> symmetry operators; sg_register_if_unknown adds to it.
sym_dict = {}
for ops, names in _tables['sym_base']:
    for sym in names:
        sym_dict[sym] = tuple(ops)

space_group_map = _tables['space_group_map']
del _tables
