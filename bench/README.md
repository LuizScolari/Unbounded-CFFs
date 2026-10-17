# Growth benchmark

Measures how long it takes to **grow one CFF into the next** (a single embedding
step) for the three constructions, counting **only the work that produces the
new parts** of the grown matrix.

The benchmark is organised as **six tables**, one per construction/regime. Each
table is a list of **blocks (chains)**: the first row of a block is an initial
CFF generated from scratch, and every row below it is an **embedding of the row
above**. Every row shows the full parameters of that CFF plus the time to
produce it.

| Table | Construction | Regime |
|---|---|---|
| `rs-embedding` | RS (polynomial) | full `q × q` rows, `d = ⌊(q−1)/k⌋`; a step raises `k`, jumps to a larger sub-field `q`, or both |
| `rs-monotone` | RS (monotone) | `d` and `k` fixed, only `q` grows (rows `(dk+1)q`) |
| `subset-embedding` | Subset | `d` fixed per chain, `k` and `S` both grow |
| `subset-nested` | Subset | nested (Sperner) family: `d = 1`, `k = ⌊S/2⌋` (maximum columns `C(S,k)`) |
| `subset-monotone` | Subset | monotone family: `d` and `k` fixed, only `S` grows |
| `crs-embedding` | CRS | one prime power added per step, `n` always maximal |

The benchmark never touches the filesystem: every step is described purely by
parameters, so **no file read/write time is ever included**. It also excludes
the re-computation of the old submatrix wherever the construction lets us
separate it. The timed region uses the construction's real inner loop (the
library functions are reused, not reimplemented), and timing is wall-clock
(`CLOCK_MONOTONIC`) over the real multi-threaded execution — i.e. the actual
latency to grow the family.

## Time columns

**RS reports two times**, mirroring the two ways the matrix gets filled:

| Column | What it measures |
|---|---|
| `Inv-CFF(s)` | The whole step: setup (sub-field partition, polynomials, `(x,y)` combinations) **+ building the inverted evaluation index** + filling the new blocks + concatenating them into the grown matrix. |
| `CFF(s)` | The same growth **with the inverted evaluation index already computed**: filling + concatenation only. |

The difference between the two columns is exactly what building the index
costs. The split comes from the library itself
(`generate_new_cff_blocks_timed`), so both numbers describe the real code path.

**Subset and CRS fill the matrix directly** — there is no index to precompute —
so they have a single `CFF(s)` column with the time to produce the new parts.

## What "new parts" means, per construction

| Construction | New parts that are timed | Excluded |
|---|---|---|
| **RS** | Generation of the three new blocks (`old_new`, `new_old`, `new_new`) via `generate_new_cff_blocks_timed` **+** concatenation (placing only those new blocks into the final matrix). | Reading the old matrix, copying the old matrix into the result, writing the output, allocating the destination buffer. |
| **CRS** | Setting the 1s **only in the new positions**: for the old prime powers, the new columns `[n_old, n_new)`; for the new prime powers, all columns. Reuses `generate_block` with a per-block `col_start`. | Re-setting the 1s of the old top-left block, allocation (which zero-fills the whole matrix), file I/O. |
| **Subset** | Generating **only the new columns** `[n_old, n_new)` with the same inner loop as `generate_cff`. | Old columns (hence old rows), allocating the destination rows, file I/O. |

For the **first row of a block** (an initial CFF) the "new parts" are the whole
matrix, so its time is the from-scratch generation.

### Why these definitions

- **RS** — `embed_cff` reads the old matrix, generates the new blocks, then
  concatenates old + new into the grown matrix. All the finite-field cost lives
  in `generate_new_cff_blocks`; the concatenation of the new blocks is the extra
  cost of assembling them. That sum is exactly "generation + concatenation of the
  new parts". The old/new block geometry is derived from the generated blocks
  themselves (old rows = rows of `old_new`, old cols = cols of `new_old`),
  matching `embed_cff` exactly, so nothing is read from disk — and because that
  geometry is the library's own, the harness cross-checks it against the `(t,n)`
  written in the table and warns on any mismatch.

- **CRS** — the embedding does not read the old matrix; it regenerates
  everything and re-sets the 1s of the old block too. We time only the 1s of the
  new positions, which is what actually grows.

- **Subset** — the embedding reorders columns so the predecessor occupies the
  upper-left block; the new columns are therefore exactly the indices
  `[n_old, n_new)` (verified empirically: the first `t_old × n_old` block of an
  embedded matrix is bit-identical to the predecessor). Because generation is
  column-parallel, a column is the natural atomic unit, and the growth cost is
  dominated by materialising these new columns. This works for **all three**
  Subset regimes:
  - *monotone* (`k' = k`): old columns hold only old cells, so the new columns
    are exactly the new cells;
  - *nested / embedding* (`k` grows): old columns are deliberately **not**
    processed (the convention is that old columns/rows are excluded), and each
    new column is materialised in full.

## Parameters

The chains are plain arrays in the `build_*` functions of
[`bench_plan.c`](bench_plan.c) — that file is the single place to edit to
benchmark other families.

- **`rs-embedding`** — three chains, one per characteristic, each staying inside
  its field tower so every `q` is a sub-field of the next:
  `(q,k) = (2,1) (4,1) (4,2) (16,2) (16,3) (16,4)`,
  `(3,1) (3,2) (9,2) (9,3) (9,4) (9,5)` and `(5,1) (5,2) (25,2) (25,3)`.
  The strength is always `d = ⌊(q−1)/k⌋`, the largest the full `q × q` block
  supports (which is also the condition under which the library uses it).
- **`rs-monotone`** — `d,k` fixed, `q` grows: `d=1,k=1: 2,4,16,256`;
  `d=2,k=1: 3,9,81`; `d=1,k=2: 3,27`; `d=2,k=2: 5,25`.
- **`subset-embedding`** — `d=1: (k,S) = (5,10) (8,16) (11,22) (13,26)`;
  `d=2: (4,10) (7,16) (10,22)`; `d=3: (5,12) (8,18) (10,22)`.
- **`subset-nested`** — `d=1`, `k=⌊S/2⌋`: `S = 12,16,20,24,28` and
  `S = 11,15,19,23,27`.
- **`subset-monotone`** — `d=1,k=10: S=20,24,28,32`; `d=2,k=8: S=16,20,24,28`;
  `d=3,k=6: S=12,16,20,24`.
- **`crs-embedding`** — `d=1` from `{2,3}` adding `5,7,11,13,17,19`; `d=2` from
  `{2,3,5,7}` adding `11,13,17,19,23,29,31,37`; `d=3` from `{2,3,5,7,11}` adding
  `13,17,19,23,29,31`. `n` is always the maximum allowed,
  `n = ⌊(∏q)^{1/d}⌋`.

## Methodology

- **Timer**: `clock_gettime(CLOCK_MONOTONIC)`, wall-clock seconds.
- **Repetitions**: `--reps N` (default **100**) per row. Reported: **mean**; the
  CSV also carries min, max, sample std-dev and median.
- **Time budget**: `--budget S` (default **20 s**) stops repeating a row once it
  has spent `S` seconds, so a row costing minutes still gets one measurement
  while cheap rows get the full 100. The actual repetition count is reported per
  row (`reps` column of the CSV and of the console tables).
- **Warm-up**: `--warmup W` (default **5**) discarded runs per row, capped at
  `budget/10` seconds so heavy rows only pay one. This matters: the first
  parallel region of a process costs milliseconds (FLINT/GLib context init,
  OpenMP/pthread pool spin-up, first-touch page faults), which would otherwise
  swamp the sub-millisecond rows.
- **Memory cap**: `--max-bytes B` skips any row whose **estimated peak memory**
  would exceed `B` (default `4G`), a guard against OOM. Skipped rows appear as
  `—` in the LaTeX tables. For Subset and CRS the estimate is the bit-packed
  matrix; for RS the matrix is the small part — the inverted evaluation index
  holds one polynomial index per (evaluation point, polynomial) pair inside
  GLib containers, measured at roughly 70 bytes per entry, so the estimate is
  dominated by `70·q·n` (the largest RS row of the default plan, `q=16, k=4`,
  peaks around 1.2 GB).
- **Noise**: rows in the microsecond range are dominated by parallel-region
  overhead and are sensitive to whatever else the machine is doing. Run the
  benchmark on an otherwise idle machine, and consider
  `gen_tables.py --stat median` for tables where those rows matter.

## Build & run

```bash
make bench                 # produces ./cff_bench (does not rebuild ./cff)

./cff_bench                        # all six tables, up to 100 reps per row
./cff_bench --list                 # list the plan without running it
./cff_bench --only rs              # rs-embedding + rs-monotone
./cff_bench --only subset-nested   # a single table
./cff_bench --reps 5 --budget 3    # quick pass
./cff_bench --csv bench/results.csv
./cff_bench --max-bytes 300M       # skip the heaviest rows
```

`--only` accepts a table id (`rs-embedding`, `rs-monotone`, `subset-embedding`,
`subset-nested`, `subset-monotone`, `crs-embedding`), a construction group
(`rs`, `subset`, `crs`) or `all`.

## Output

The console prints the six tables (parameters + times), and `--csv` writes one
row per plan row with everything:

```
table, construction, block, step, q, k, d, S, t, n, t_prev, n_prev, reps, ran,
inv_mean_s, inv_min_s, inv_std_s, inv_median_s,
cff_mean_s, cff_min_s, cff_std_s, cff_median_s,
concat_mean_s, new_cells, matrix_MB
```

`inv_*` is empty for Subset and CRS (no inverted index); `q` holds the field
size for RS and the prime power(s) for CRS, and is quoted because the CRS
initial row lists several (`"2,3,5,7"`). `ran=0` marks a row skipped by the
memory cap.

## LaTeX tables

`gen_tables.py` turns the CSV into the six LaTeX tables, each with its blocks
separated by `\addlinespace`:

```bash
# one shot: benchmark -> CSV -> bench/tables.tex
make tables
make tables REPS=5                 # quicker pass

# or by hand:
./cff_bench --reps 100 --csv bench/results.csv
python3 bench/gen_tables.py bench/results.csv --out bench/tables.tex
python3 bench/gen_tables.py bench/results.csv --standalone -o doc.tex  # compilable
python3 bench/gen_tables.py bench/results.csv --stat median            # robust stat
python3 bench/gen_tables.py bench/results.csv --only rs-embedding      # one table
python3 bench/gen_tables.py bench/results.csv --decimals 4
```

The bare output (no `--standalone`) is meant to be `\input`; it needs
`\usepackage{booktabs}` and `\usepackage{amsmath}` in your preamble (emitted as
a comment at the top of the file, and in full with `--standalone`).

Because the tables are driven entirely by the CSV, re-running the benchmark with
different repetitions or parameter chains regenerates them automatically.

## Relationship to the library

The benchmark reuses the library's real inner loops. To make them reachable from
`bench/`, a few functions were given external linkage (no behaviour change to
`./cff`):

- RS: `generate_new_cff_blocks`, `generate_new_cff_blocks_timed`, `free_matrix`
  (declared in `src/RS/cff_builder.h`). The `_timed` variant is
  `generate_new_cff_blocks` with an optional `rs_phase_times` out-parameter
  reporting the setup / inverted-index / fill split; the original is now a
  one-line wrapper around it, so `./cff` runs exactly the same code.
- CRS: `generate_block` (declared in `src/CRS/generate.h`).
- Subset: `chain_col`, `next_comb_init`, `next_comb` (declared in
  `src/Subset/generate_cff.h`).
