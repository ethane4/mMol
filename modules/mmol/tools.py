"""Python side of mMol's analysis tabs: Ramachandran plot, Align, Mutagenesis.

Each call that the Swift panel reads back writes JSON to `path` (feedback
lines are too short for the payloads).
"""

import json

from pymol import cmd


def _write(path, data):
    with open(path, 'w') as f:
        json.dump(data, f)


# --- Ramachandran -----------------------------------------------------------


def rama(obj, path):
    """phi/psi of every protein residue of `obj` (current state)."""
    try:
        angles = cmd.phi_psi('(%s) and polymer.protein and name CA' % obj)
        info = {}
        cmd.iterate('(%s) and polymer.protein and name CA' % obj,
                    'info[(model, index)] = (chain, resi, resn)',
                    space={'info': info})
        res = [
            dict(chain=info[k][0],
                 resi=info[k][1],
                 resn=info[k][2],
                 phi=v[0],
                 psi=v[1]) for k, v in angles.items() if k in info
        ]
        res.sort(key=lambda r: (
            r['chain'],
            int(''.join(c for c in r['resi'] if c.isdigit() or c == '-') or 0)))
        _write(path, {'res': res})
    except Exception as e:
        _write(path, {'error': str(e)})


# --- Align --------------------------------------------------------------------


def align(mobile, mobile_region, target, target_region, method, ca_only, path):
    """Superpose `mobile` onto `target` using only the two regions, then report
    the region RMSD (atoms used for the fit) and the whole-object Cα RMSD."""

    def sel(obj, region):
        s = '(%s) and (%s)' % (obj, region.strip() or 'all')
        return s + ' and name CA' if ca_only else s

    try:
        m, t = sel(mobile, mobile_region), sel(target, target_region)
        for s, what in ((m, 'mobile'), (t, 'target')):
            if not cmd.count_atoms(s):
                raise ValueError('No %satoms in the %s region.' %
                                 ('Cα ' if ca_only else '', what))
        out = {'method': method}
        if method == 'cealign':
            r = cmd.cealign(t, m)
            out.update(rmsd=r['RMSD'], atoms=r['alignment_length'])
        else:
            r = getattr(cmd, method)(m, t)
            out.update(rmsd=r[0],
                       atoms=r[1],
                       rmsd_before=r[3],
                       atoms_before=r[4])
        # Whole objects, sequence-matched Cα pairs, measured without moving.
        tot = cmd.align('(%s) and name CA' % mobile,
                        '(%s) and name CA' % target,
                        cycles=0,
                        transform=0)
        out.update(total_rmsd=tot[0], total_atoms=tot[1])
        print(
            ' Align: region RMSD %.3f A over %d atoms; whole-object CA RMSD %.3f A over %d pairs.'
            %
            (out['rmsd'], out['atoms'], out['total_rmsd'], out['total_atoms']))
    except Exception as e:
        out = {'error': str(e) or type(e).__name__}
    _write(path, out)


# --- Mutagenesis ------------------------------------------------------------


def _wizard():
    from pymol.wizard.mutagenesis import Mutagenesis
    w = cmd.get_wizard()
    return w if isinstance(w, Mutagenesis) else None


def mut_status(path, error=None):
    w = _wizard()
    out = {'error': error}
    if w is not None and w.status == 1:
        n = cmd.count_states('mutation') if 'mutation' in cmd.get_names() else 0
        out.update(residue=w.res_text,
                   mode=w.mode,
                   states=n,
                   scores=list(w.bump_scores),
                   state=cmd.get_state())
    _write(path, out)


def mut_select(selection, path):
    """Start (or move) the wizard onto the first residue of `selection`."""
    from pymol.wizard.mutagenesis import Mutagenesis
    try:
        if not cmd.count_atoms(selection):
            return mut_status(
                path, 'Select a residue first (click one in the viewer).')
        w = _wizard()
        if w is None:
            w = Mutagenesis(_self=cmd)
            w.rep = 'sticks'
            cmd.set_wizard(w)
        cmd.select('_mmol_mut',
                   'byres first (%s)' % selection)  # do_select deletes it
        w.do_select('_mmol_mut')
        mut_status(path)
    except Exception as e:
        mut_status(path, str(e) or type(e).__name__)


def mut_mode(mode, path):
    w = _wizard()
    if w is not None:
        w.set_mode(mode)
    mut_status(path)


def mut_state(state, path):
    w = _wizard()
    if w is not None:
        cmd.frame(state)
        w.do_state(state)
    mut_status(path)


def mut_apply(path):
    w = _wizard()
    if w is not None:
        w.apply()
    mut_status(path)


def mut_end():
    w = _wizard()
    if w is not None:
        w.clear()
        cmd.set_wizard()
