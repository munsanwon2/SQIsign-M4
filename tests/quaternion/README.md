# Quaternion arithmetic regressions

Run these tests with 64-bit Linux or WSL, GCC or Clang and a GNU-compatible
linker. The runner uses the checked 34/53/69-limb radix-32 exports and the
compiler's UBSan runtime.

```sh
bash tests/quaternion/run.sh
CC=clang bash tests/quaternion/run.sh 3
```

- `normeq_regression.c` covers canonical residues and signed overflow in
  the congruence checks ported from compact-SQIsign PR #3. The test fixes
  candidate selection and checks the exact response norm and unpublished
  output on failure.
- `ml2_regression.c` checks dependent rows, zero generators, rank and
  exact lattice equality. It tests the rational Lovász and size-reduction
  conditions at δ = 43/50 and η = 11/20. It also compares on-demand Gram
  recomputation with an independent exact snapshot update.
- `compact_regression.c` selects the portable 32-bit-half integer path.
  It checks signed capacity, both digit widths, DPE halfway rounding and
  large exponents, ML2's 512-bound and rounding-mode restoration, and
  the direct sampler's small-ball support and failure statuses. It also
  checks guard reentry, workspace exhaustion, leak reporting and recovery.

The host tests build into a temporary directory. End-to-end signing and
target checks are described in [VALIDATION.md](../../VALIDATION.md).
