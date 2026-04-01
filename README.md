# Unbounded CFFs Implementation

This project implements the generation and embedding (expansion) of **Cover-Free Families (CFFs)** using constructions based on polynomials over finite fields. The code is written in C and optimized for performance.

## What are CFFs?

**Cover-Free Families (CFFs)** are combinatorial structures which can be represented as binary matrices. This matrix is considered a $d$-CFF if the union of any $d$ columns does not completely cover any other remaining column.

They have important applications in:

  * Group Testing
  * Cryptography and key distribution
  * Digital signatures

## Implemented Constructions

The project focuses on algebraic construction using polynomials over finite fields ($\mathbb{F}_q$).

### Construction of Polynomials over Finite Fields - Embedding CFFs

In this construction, matrix rows are indexed by pairs of elements $(x, y) \in B \times \mathbb{F}_q$ where $B \in \mathbb{F}_q$ , and columns are indexed by polynomials of degree up to $k$. A position in the matrix is set to $1$ if the corresponding polynomial evaluated at $x$ yields $y$ (i.e., $P(x) = y$), and 0 otherwise.

### Monotone CFFs

The Monotone construction is a variation designed for embedding operations. During expansion, rows are indexed by pairs $(x, y)$ where $x$ is restricted to a fixed subset $B \subseteq \mathbb{F}_q$ corresponding to the smaller base field. The parameters $d$ and $k$ are always kept constant throughout the embedding process—only the field size $q$ changes.

## Hash Table Optimizations

Naive generation of polynomial CFFs can be computationally expensive due to the need to evaluate many polynomials at many points.

To optimize this process, this project uses **Hash Tables (GHashTable from GLib)** to create an inverted index of evaluations:

1.  Instead of repeatedly recalculating $P(x)$, we precompute and store which polynomials yield a value $y$ for a given $x$.
2.  This transforms the search into an average constant-time $O(1)$ operation for filling the matrix, drastically speeding up the generation of large instances.

## Prerequisites

To compile this project, you will need the following libraries installed on your system:

  * **GCC** or **Clang**
  * **Make**
  * **FLINT 3.3.1** (Fast Library for Number Theory)
  * **GMP** (GNU Multiple Precision Arithmetic Library)
  * **GLib 2.0**
  * **OpenMP** (libomp)

### Supported Platforms

This project supports **Linux** and **macOS**. The Makefile automatically detects the operating system and configures the appropriate compiler flags.

| Platform | Compiler | Notes |
|----------|----------|-------|
| Linux | GCC | Standard configuration |
| macOS | Clang | Requires Homebrew for dependencies |

#### Installing Dependencies

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install build-essential libgmp-dev libglib2.0-dev libomp-dev
# FLINT 3.3.1 must be compiled from source (see flintlib.org)
```

**macOS (Homebrew):**
```bash
brew install gmp glib libomp flint
```

## How to Run

### 1\. Compilation

Use the provided `Makefile` to compile the project:

```bash
make
```

### 2\. Execution

The `generate_cff` executable supports different operation modes. Arguments vary according to the construction type (`p` for initial and embedding CFFs, `m` for monotone CFFs) and action (`f` to generate initial CFF, `g` for embedding/expansion).

**General Syntax:**
`./generate_cff <type> <action> [parameters...]`

**Block Size Option (Initial or Embedding CFFs only):**
  * `m` - Minimum block size, calculated as: $(dk+1)q$
  * `f` - Full block size, calculated as: $q^2$
  
Note: When $d = (q-1)/k$, the minimum block size $(dk+1)q$ can yield fewer rows than the full size $q^2$. Example: $\mathbb{F}_4$, $k=2$ → $12 < 16$ rows. If $d < (q-1)/k$, $(dk+1)q$ is always considered.

#### Examples:

  * **Generate an initial CFF from scratch (`p f`):**

      * Parameters: `<m|f>` (block size), `d` (target $d$), `q` (target $q$), `k` (target $k$).

    <!-- end list -->

    ```bash
    ./generate_cff p f m 2 3 1
    ```

  * **Embedding CFF (Expansion) (`p g`):**

      * Expands an existing CFF to a larger field.
      * Parameters: `<m|f>` (block size), `cff_file`, `d` (target $d$), `q` (target $q$), `k` (target $k$).

    <!-- end list -->

    ```bash
    ./generate_cff p g m CFFs/cff_input.txt 2 9 1
    ```

  * **Monotone CFF (Expansion) (`m g`):**

      * Parameters: `cff_file`, `d` (fixed $d$), `q` (target $q$), `k` (fixed $k$).

    <!-- end list -->

    ```bash
    ./generate_cff m g CFFs/cff_input.txt 2 9 1
    ```

Output files will be generated in the `CFFs/` folder.

## Benchmark

The project includes an automated benchmark system to measure the execution time of CFF generation. The benchmarks measure two metrics:

- **Time 1**: Total time for inverted index + CFF generation (+ concatenation for embedding)
- **Time 2**: Time for CFF matrix generation + concatenation only

### Quick Start

```bash
cd benchmark
make
./run_benchmarks.sh
```

For more details on how to run the benchmarks, configure iterations, and interpret results, see the [benchmark/README.md](benchmark/README.md).

-----

## Supervision

This research is being supervised by **Professor Thaís Bardini Idalino** ([Website](https://thaisidalino.github.io), [Google Scholar](https://scholar.google.com/citations?user=hS_R7ZsAAAAJ&hl)).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
