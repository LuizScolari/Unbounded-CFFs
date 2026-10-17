#!/usr/bin/env python3
"""
gen_tables.py - Generate the LaTeX growth-benchmark tables from a cff_bench CSV.

One table per construction/regime (six in total):

    rs-embedding      RS polynomial, full q x q rows, d = floor((q-1)/k)
    rs-monotone       RS monotone: d, k fixed, only q grows
    subset-embedding  Subset: d fixed, k and S grow
    subset-nested     Subset nested (Sperner): d = 1, k = floor(S/2)
    subset-monotone   Subset monotone: d, k fixed, S grows
    crs-embedding     CRS: one prime power added per step, n maximal

Every table lists the parameters of each CFF followed by the growth time. Rows
are grouped in blocks (chains): the first row of a block is an initial CFF
generated from scratch and each row below it is an embedding of the row above.
RS carries two time columns -- Inv-CFF(s), which includes building the inverted
evaluation index, and CFF(s), the same growth with the index already computed;
the other constructions fill the matrix directly and have a single CFF(s).

Usage:
    ./cff_bench --csv bench/results.csv                     # produce the data
    python3 bench/gen_tables.py bench/results.csv           # -> stdout
    python3 bench/gen_tables.py bench/results.csv --out bench/tables.tex
    python3 bench/gen_tables.py bench/results.csv --standalone -o doc.tex
    python3 bench/gen_tables.py bench/results.csv --stat median
    python3 bench/gen_tables.py bench/results.csv --only rs-embedding

The bare output needs \\usepackage{booktabs} in your preamble (emitted as a
comment at the top of the file, and in full with --standalone).
"""

import argparse
import csv
import sys

MISSING = r'\textemdash'


# --------------------------------------------------------------- table specs

# id -> (title, parameter columns, time columns, caption)
# A parameter column is (CSV field, LaTeX header, alignment).
PARAM_RS = [('q', '$q$', 'c'), ('k', '$k$', 'c'), ('d', '$d$', 'c'),
            ('t', '$t$', 'r'), ('n', '$n$', 'r')]
PARAM_SUB = [('S', '$S$', 'c'), ('k', '$k$', 'c'), ('d', '$d$', 'c'),
             ('t', '$t$', 'r'), ('n', '$n$', 'r')]
PARAM_CRS = [('q', '$q$', 'c'), ('d', '$d$', 'c'),
             ('t', '$t$', 'r'), ('n', '$n$', 'r')]

TIME_RS = [('inv', 'Inv-CFF(s)'), ('cff', 'CFF(s)')]
TIME_ONE = [('cff', 'CFF(s)')]

TABLES = [
    ('rs-embedding', PARAM_RS, TIME_RS,
     r'Embedding family from the polynomial (Reed--Solomon) construction.'),
    ('rs-monotone', PARAM_RS, TIME_RS,
     r'Monotone family from the polynomial (Reed--Solomon) construction.'),
    ('subset-embedding', PARAM_SUB, TIME_ONE,
     r'Embedding family from the subset construction.'),
    ('subset-nested', PARAM_SUB, TIME_ONE,
     r'Nested family from the subset construction.'),
    ('subset-monotone', PARAM_SUB, TIME_ONE,
     r'Monotone family from the subset construction.'),
    ('crs-embedding', PARAM_CRS, TIME_ONE,
     r'Embedding family from the Chinese Remainder Sieve (CRS) construction.'),
]


# --------------------------------------------------------------- formatting

def group(value):
    """'1048576' -> '1{,}048{,}576' (thousands separators, LaTeX-safe)."""
    s = str(value).strip()
    if not s.lstrip('-').isdigit():
        return s.replace(',', '{,}')          # e.g. the CRS prime-power list
    neg, s = s.startswith('-'), s.lstrip('-')
    out = []
    for i, ch in enumerate(s):
        if i and (len(s) - i) % 3 == 0:
            out.append('{,}')
        out.append(ch)
    return ('-' if neg else '') + ''.join(out)


def fmt_time(row, prefix, stat, decimals):
    value = row.get(f'{prefix}_{stat}_s', '')
    if row.get('ran') != '1' or value in (None, ''):
        return MISSING
    return f'{float(value):.{decimals}f}'


# ------------------------------------------------------------- table emitter

def emit_table(table_id, params, times, caption, rows, stat, decimals):
    rows = [r for r in rows if r['table'] == table_id]
    if not rows:
        return ''

    align = ''.join(a for _, _, a in params) + 'r' * len(times)
    header = [h for _, h, _ in params] + [r'\textbf{%s}' % h for _, h in times]

    lines = [r'\begin{table}[htbp]', r'  \centering',
             r'  \begin{tabular}{%s}' % align, r'  \toprule',
             '  ' + ' & '.join(header) + r' \\', r'  \midrule']

    prev_block = None
    for r in rows:
        if prev_block is not None and r['block'] != prev_block:
            lines.append(r'  \addlinespace')
        prev_block = r['block']

        cells = [group(r[field]) for field, _, _ in params]
        cells += [fmt_time(r, prefix, stat, decimals) for prefix, _ in times]
        lines.append('  ' + ' & '.join(cells) + r' \\')

    reps = max((int(r['reps']) for r in rows if r.get('reps')), default=0)
    lines += [r'  \bottomrule', r'  \end{tabular}',
              r'  \caption{%s}' % caption.replace('@reps@', str(reps))
                                        .replace('@stat@', stat),
              r'  \label{tab:%s}' % table_id,
              r'\end{table}']
    return '\n'.join(lines)


PREAMBLE_NOTE = r"""% Generated by bench/gen_tables.py -- CFF growth benchmark tables.
% Requires, in your document preamble:
%   \usepackage{booktabs}
"""

STANDALONE_HEAD = r"""\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage[margin=2cm]{geometry}
\begin{document}
"""

STANDALONE_TAIL = r'\end{document}'


def main():
    ap = argparse.ArgumentParser(
        description='Generate the LaTeX growth tables from a cff_bench CSV.')
    ap.add_argument('csv', help='CSV produced by ./cff_bench --csv')
    ap.add_argument('-o', '--out', help='output .tex (default: stdout)')
    ap.add_argument('--standalone', action='store_true',
                    help='emit a full compilable document instead of bare tables')
    ap.add_argument('--stat', default='mean', choices=['mean', 'median', 'min'],
                    help='statistic shown in the time columns (default mean; '
                         'median is the robust choice for sub-millisecond rows)')
    ap.add_argument('--decimals', type=int, default=6,
                    help='decimal places in the time columns (default 6)')
    ap.add_argument('--only', help='emit a single table by id (e.g. rs-embedding)')
    args = ap.parse_args()

    with open(args.csv, newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        sys.exit('Empty CSV.')

    specs = TABLES if not args.only else [t for t in TABLES if t[0] == args.only]
    if not specs:
        sys.exit(f'Unknown table id: {args.only}')

    tables = [emit_table(tid, params, times, caption, rows, args.stat, args.decimals)
              for tid, params, times, caption in specs]
    tables = [t for t in tables if t]
    if not tables:
        sys.exit('No measured rows in the CSV.')

    parts = [STANDALONE_HEAD if args.standalone else PREAMBLE_NOTE,
             '\n\n'.join(tables)]
    if args.standalone:
        parts.append('\n' + STANDALONE_TAIL)
    text = '\n'.join(parts) + '\n'

    if args.out:
        with open(args.out, 'w') as f:
            f.write(text)
        print(f'Wrote {args.out} ({len(tables)} table(s)).', file=sys.stderr)
    else:
        sys.stdout.write(text)


if __name__ == '__main__':
    main()
