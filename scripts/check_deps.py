#!/usr/bin/env python3
"""Check that the C++ components under src/ only depend downward.

Each folder of ``src/pymol`` (and ``src/ov``) is a component. The components
are ordered in layers; a component may include headers from its own folder and
from lower layers only. Components in the same layer may not include each
other, except in the core layer, whose components still depend on each other
and are untangled one at a time (see CLAUDE.md). A component's
``*_internal.h`` headers are private to it. Third-party and system headers
are not checked.

Usage::

    scripts/check_deps.py            # report violations, exit 1 if any
    scripts/check_deps.py --graph    # also print each component's uses
"""

import argparse
import collections
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'src')

# Lowest layer first.
LAYERS = [
    ['ov'],
    ['platform'],
    ['base'],
    ['pyconv'],
    ['settings'],
    ['geometry'],
    ['align', 'image', 'symmetry'],
    ['color', 'fields', 'glyph'],
]
CORE = [
    'app', 'chem', 'executive', 'graphics', 'io', 'objects', 'python', 'ray',
    'reps', 'scene', 'select', 'text', 'ui'
]

INCLUDE = re.compile(r'^\s*#\s*include\s+"([^"]+)"', re.M)


def component_of(path):
    """Return the component a path under ``src/`` belongs to.

    Parameters
    ----------
    path : str
        Path relative to ``src/``, such as ``pymol/scene/scene.h``.

    Returns
    -------
    str or None
        The component name (``scene``, ``ov``), or None for files outside any
        component (the shared ``pymol/*.h`` utility headers).
    """
    parts = path.split('/')
    if parts[0] == 'ov':
        return 'ov'
    if parts[0] == 'pymol' and len(parts) > 2:
        return parts[1]
    return None


def layer_index():
    """Map each component to its layer number.

    Returns
    -------
    dict of str to int
        Layer of each component, 0 for the lowest; the core is the top layer.
    """
    index = {c: i for i, layer in enumerate(LAYERS) for c in layer}
    index.update({c: len(LAYERS) for c in CORE})
    return index


def scan():
    """Collect the cross-component includes under ``src/``.

    Returns
    -------
    dict of (str, str) to list of str
        For each (component, used component) pair, the ``file:line`` sites of
        the includes.
    """
    edges = collections.defaultdict(list)
    for root, _, files in os.walk(SRC):
        for name in files:
            if not name.endswith(('.h', '.cc', '.mm')):
                continue
            path = os.path.join(root, name)
            rel = os.path.relpath(path, SRC)
            src = component_of(rel)
            if src is None:
                continue
            with open(path, encoding='utf-8', errors='replace') as f:
                text = f.read()
            for m in INCLUDE.finditer(text):
                dst = component_of(m.group(1))
                if dst and dst != src:
                    line = text.count('\n', 0, m.start()) + 1
                    site = '%s:%d' % (rel, line)
                    if m.group(1).endswith('_internal.h'):
                        site += ' (private header %s)' % m.group(1)
                    edges[(src, dst)].append(site)
    return edges


def violations(edges, index):
    """Find the includes that break the layering.

    Parameters
    ----------
    edges : dict of (str, str) to list of str
        Cross-component includes, as `scan` returns them.
    index : dict of str to int
        Layer of each component, as `layer_index` returns it.

    Returns
    -------
    list of str
        One ``file:line: from -> to (why)`` message per offending include.
    """
    out = []
    for (src, dst), sites in sorted(edges.items()):
        if src not in index or dst not in index:
            out += [
                '%s: unknown component %r' %
                (s, src if src not in index else dst) for s in sites
            ]
            continue
        upward = index[dst] > index[src]
        sibling = index[dst] == index[src] and src not in CORE
        if upward or sibling:
            why = 'higher layer' if upward else 'same layer'
            out += ['%s: %s -> %s (%s)' % (s, src, dst, why) for s in sites]
        out += [
            '%s: %s -> %s' % (s, src, dst)
            for s in sites
            if 'private header' in s and not (upward or sibling)
        ]
    return out


def main():
    """Run the check; return the process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--graph',
                        action='store_true',
                        help="print each component's dependencies")
    args = parser.parse_args()
    index = layer_index()
    edges = scan()
    if args.graph:
        uses = collections.defaultdict(set)
        for src, dst in edges:
            uses[src].add(dst)
        for c in sorted(uses, key=lambda c: (index.get(c, 99), c)):
            print('%-10s %s' % (c, ' '.join(sorted(uses[c]))))
    bad = violations(edges, index)
    for msg in bad:
        print(msg)
    print('check_deps: %d violation(s)' % len(bad), file=sys.stderr)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
