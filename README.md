# Unbounded Cover-Free Families

A C library gathering three constructions of **unbounded Cover-Free Families (CFFs)**
behind a single command-line tool. Each construction can generate an initial CFF
and **embed** an existing one into a larger family, keeping the predecessor as a
submatrix (which is what makes the families *unbounded*).

## What is a CFF?

A **d-CFF(t, n)** is a *t × n* binary matrix such that, for any *d + 1* columns,
there is a row with a `1` in one of them and `0` in all the others. Equivalently,
the union of any *d* columns never covers another column. The ratio *n/t* (or its
growth) measures how efficient the family is.

CFFs are used in group testing, cryptography and key distribution, digital
signatures, and combinatorics.

## Repository layout

```
src/
  main.c        dispatcher: first argument selects the construction
  RS/           Construction of polynomials over finite fields    
  CRS/          Chinese Remainder Sieve construction           
  Subset/       Subset construction            
test/           CFF verifiers: verify.c (d-CFF property) and verify.py (row/col sums)
Makefile        single build → ./cff
CFFs/           output matrices (created at runtime, git-ignored)
```

## Build

A single `make` produces one binary, `./cff`:

```bash
make
make clean
```

The tool works on **Linux** (GCC) and **macOS** (Clang + Homebrew); the Makefile
detects the OS automatically.

### Dependencies

Because the binary links all three constructions together, all dependencies are
required to build it:

| Library | Used by | Install (Debian/Ubuntu) | Install (macOS) |
|---------|---------|-------------------------|-----------------|
| GCC / Clang + Make | all | `build-essential` | Xcode CLT |
| pkg-config | build | `pkg-config` | `brew install pkg-config` |
| GLib 2.0 | RS | `libglib2.0-dev` | `brew install glib` |
| GMP | RS | `libgmp-dev` | `brew install gmp` |
| FLINT 3.3.1 | RS | build from source ([flintlib.org](https://flintlib.org)) | `brew install flint` |
| OpenMP | Subset (RS) | `libomp-dev` | `brew install libomp` |
| pthreads | CRS | included | included |

## Usage

The **first argument is always the construction** (`RS`, `Subset` or `CRS`,
case-insensitive); the remaining arguments belong to that construction:

```bash
./cff <construction> <args...>
```

Run `./cff` with no arguments for the top-level help, or `./cff <construction>`
for a construction's own help. Every construction writes its result to `CFFs/`,
named automatically as `d-CFF(t,n).txt`.

---

### RS (Reed-Solomon) - construction of polynomials over finite fields

Rows are indexed by pairs *(x, y) ∈ B × 𝔽_q* and columns by polynomials of degree
up to *k* over the finite field 𝔽_q; entry *(x, y), P* is `1` iff *P(x) = y*.
Evaluations are indexed with a GLib hash table, turning matrix filling into an
average *O(1)* lookup. The **monotone** variant (`m`) keeps *d* and *k* fixed and
only grows *q* during embedding.

```bash
# Generate an initial CFF from scratch:    RS p f <m|f> <d> <q> <k>
./cff RS p f m 2 3 1

# Embed (expand) an existing CFF:          RS p g <m|f> <cff_file> <d> <q> <k>
./cff RS p g f "CFFs/2-CFF(9,9).txt" 4 9 2

# Monotone embedding:                      RS m g <cff_file> <d> <q> <k>
./cff RS m g "CFFs/2-CFF(9,9).txt" 2 9 1
```

Parameters:

- `p` / `m` — construction type: polynomial or monotone.
- `f` / `g` — action: `f` generate an initial CFF, `g` embed/expand (monotone only supports `g`).
- `<m|f>` — block (row) size, for the polynomial type only:
  - `m` — minimum, `(dk+1)·q` rows.
  - `f` — full, `q²` rows.
  - When *d = (q−1)/k*, the minimum can be smaller than the full size (e.g. 𝔽₄, *k = 2* → 12 < 16 rows).
- `d`, `q`, `k` — target strength, field size (a prime power) and polynomial degree.

### Subset construction

Builds **M(d, k, S)** with *t = C(S, d)* rows (one per *d*-subset of {0,…,S−1}),
*n = C(S, k)* columns (one per *k*-subset), and entry `1` iff the *d*-subset is
contained in the *k*-subset. Embedding into M(d, k′, S′) with *k ≤ k′ < S′* and
*S′ > S* reorders columns (colex ranks) so every predecessor lands in the
upper-left block. Columns are processed in 64-wide bit-packed blocks, one per
OpenMP thread.

```bash
# Generate M(d,k,S):                               Subset generate <d> <k> <S>
./cff Subset generate 1 5 10                       # → CFFs/1-CFF(10,252).txt

# Embed into a larger one:                         Subset embed <old_file> <k'> <S'>
./cff Subset embed "CFFs/1-CFF(10,252).txt" 7 14   # → CFFs/1-CFF(14,3432).txt
```

Parameters:

- `generate` requires `0 < d < k < S`.
- `embed` requires `k ≤ k'`, `S < S'`, and `k' < S'`.
- keeping `k` fixed across embeddings (`k' = k`) yields a monotone family.
- `d` — strength; keep `d = 1` and `k = ⌊S/2⌋` at each step for the optimal
  (Sperner) nested family with ratio *ρ(n) ~ n/log n*.

### CRS — Chinese Remainder Sieve construction

Given *k* prime powers *p₁,…,p_k* with distinct primes (pairwise coprime moduli)
and *∏ pⱼ ≥ nᵈ*, each prime power contributes a block of *pⱼ* rows; entry
*((j, y), c)* is `1` iff *c ≡ y (mod pⱼ)*.
The result is a *d-CFF(t, n)* with *t = Σ pⱼ*. Embedding adds new prime powers
and/or a larger *n* while preserving the original matrix. Row blocks are generated
in parallel with pthreads over a bit-packed matrix.

```bash
# Generate:                                     CRS generate <d> <n> <q1> <q2> ... <qk>
./cff CRS generate 1 30 2 3 5                   # → CFFs/1-CFF(10,30).txt

# Embed:                                        CRS embed <old_file> <d'> <n'> <q_new1> ...
./cff CRS embed "CFFs/1-CFF(10,30).txt" 1 50 9  # → CFFs/1-CFF(19,50).txt
```

Parameters:

- `d` — target strength.
- `n` — number of items (columns).
- `q1 … qk` — the prime powers (with distinct primes), one per row block.
- For `embed`, `d'` / `n'` are the new strength / column count and the trailing
  values are the prime powers to add.

## File format

Output files live in `CFFs/` as `d-CFF(t,n).txt`. The first line holds the
parameters (the exact fields depend on the construction) and each remaining line
is one binary matrix row. The Subset and CRS embeddings also record the
embedding chain in the header so predecessors can be located in the upper-left
block.

## Supervision

This research is supervised by **Professor Thaís Bardini Idalino**
([website](https://thaisidalino.github.io),
[Google Scholar](https://scholar.google.com/citations?user=hS_R7ZsAAAAJ&hl)).

## License

Licensed under the MIT License — see [LICENSE](LICENSE).
