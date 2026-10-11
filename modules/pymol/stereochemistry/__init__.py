import sys
import os
import pymol

cmd = sys.modules["pymol.cmd"]


def assign_stereo(selection='all',
                  state=-1,
                  method='',
                  quiet=1,
                  prop='stereo',
                  _self=cmd):
    '''Assign "stereo" atom property (R/S stereochemistry).

    Requires either a Schrodinger Suite installation (SCHRODINGER
    environment variable set) or RDKit (rdkit Python module).

    Parameters
    ----------
    selection : str, optional
        atom selection {default: all}
    state : int, optional
        object state {default: -1 (current)}
    method : schrodinger or rdkit, optional
        {default: try both}

    Notes
    -----
    Command-line usage::

        assign_stereo [selection [, state [, method ]]]
    '''
    raise pymol.IncentiveOnlyException()
