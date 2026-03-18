![Build](https://github.com/johnauld/probeFuzz/actions/workflows/c-cpp.yml/badge.svg)
# probeFuzz

Finds all pairs of DNA probes within a given [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) of each other.

Given a set of 50-mer nucleotide probes, `probeFuzz` reports every pair whose sequences differ at ≤ `MAXD` positions. Probes are packed into two `uint64_t` words each, and distances are computed with XOR + popcount, making the inner loop fast enough to evaluate large probe sets efficiently. OpenMP is used to parallelize the search.

Although the algorithm is O(n^2), it can process inputs of upwards of 65K probes in under 1 second, or \~240K probes within 10 seconds (on a modern desktop CPU).

---

## Requirements

- g++ >= 13 (for C++23 support)
- OpenMP (via g++ >= 13)

On Ubuntu/Debian: `sudo apt install g++ libomp-dev`  
On macOS with Homebrew: `brew install gcc libomp`

---

## Building

```bash
# Optimized release binary (default)
make

# Debug binary with AddressSanitizer + UBSan
make debug

# Remove build artifacts
make clean
```

The release build uses `-O3 -march=native -fopenmp`. The debug binary is
written to `probeFuzz_debug`.

---

## Usage

Input is read from **stdin**, one probe per line, in comma-separated format:

```
<probe_id>,<sequence>
```

Results are written to **stdout**.

```bash
./probeFuzz [maxD] < sample_data/probe3.txt
```

maxD is the optional Hamming threshold to report. The default is 10; going much above 30 can lead to performance issues.

---

## Input format

- Each line: `id,sequence` — a probe identifier and its nucleotide sequence.
- Sequences must be at least 50 bases long; only the first 50 are considered.
- Valid bases: `A C G T` (case-insensitive).
- Malformed lines (no comma, sequence too short, invalid bases) are skipped
  with a warning on stderr.

Example:

```
aaeA,GTTGCGCTGGTGAAACAGAACTCCTTCTATGTACTGGCCTATATGGAAGAAACTAAGCTG
aaeB,CGCGCAATATCAATTAATGCAACTCTGTATCAAGCATGGCGATGGTGAAGTTGTCGATAA
```

---

## Output format

Comma-separated, one matching pair per line:

```
<hamming_distance>,<id_i>,<id_j>,<seq_i>,<seq_j>
```

Sequences in the output are always lowercase. Example:

```
3,aaeA,aaeB,gttgcgctggtgaaacagaactccttctatgtactggcctatatggaagaa,cgcgcaatatcaattaatgcaactctgtatcaagcatggcgatggtgaagt
```

Note that pair ordering is not guaranteed.

---

## Sample data

Some sample input files can be found in `sample_data/`

---

## Algorithm notes

Each 50-mer is stored as two `uint64_t` words (32 bases × 2 bits in `lo_`, 18 bases × 2 bits in `hi_`). The Hamming distance is computed by XOR-ing the two words, collapsing each 2-bit result into a single mismatch flag via bitwise operations, then summing that result with `__builtin_popcountll`. This avoids any per-base branching in the inner loop.

The O(n^2) pair search is parallelized with OpenMP using dynamic scheduling (chunk size 64) to balance load across threads. Each thread buffers its output in a local `std::ostringstream`, then flushes to stdout under a mutex at the end, avoiding lock contention in the inner loop.
