# Quaternion arithmetic regressions

These tests cover the `normeq.c` and `lll/ml2.c` fixes ported from
[compact-SQIsign PR #3](https://github.com/pqc-lab-ku/compact-SQIsign/pull/3),
merged as `8f58d40328fd2e5bce33e3d388a8227bd17e83f5`.
They are separate from the DPE zero-operand change.

Run explicitly on a Linux/WSL host with GCC or Clang, GNU-compatible linker
wrapping/section collection, and the compiler's UBSan runtime:

```sh
bash tests/quaternion/run.sh       # Level I, III, V
CC=clang bash tests/quaternion/run.sh 3
```

The runner uses each level's own source export, `RADIX_32`, and existing
27/40/52-limb capacity, with fixed-precision overflow checks. It creates
temporary host binaries and removes them on exit. It does not build firmware,
change target compile flags, or enter the backend manifests. No test is added
to the automatic installation/release path.

- `normeq_regression.c` calls the production `quat_represent_integer()` with
  deterministic candidate selection via linker wrappers. Only sampling,
  primality selection, and Cornacchia's returned roots are controlled; the
  roots satisfy the norm equation exactly. Some injected Cornacchia targets
  may be composite: these are congruence-path unit tests, not prime-generation
  or end-to-end signature tests. Cases cover negative differences, both
  signed 32-bit overflow directions in each subtraction, multi-limb inputs,
  the parity swap, and rejection by either congruence. They verify the output
  coordinates and norm, and that rejection preserves the output. Direct
  signed `ibz_mod_ui()` cases check canonical residues. UBSan is enabled.
- `ml2_regression.c` exercises production ML2 and the public `quat_mlll()`
  wrapper at 4, 5, 8, and 16 generators. The core reproducer is
  `[2e0, 2e1, 2e2, e0+e1]`, where `ei` is coordinate vector `i`. Eager
  factorization after insertion reaches a dependent middle row before an
  independent last row. A test-only divisor check rejects zero/nonfinite
  Gram--Schmidt diagonals immediately. Additional cases include leading and
  trailing zeros, dependencies, all-zero input, and ordinary full-rank input.
  Checks cover rank, exact integer-lattice equality (membership plus index),
  size reduction/Lovasz inequalities, input preservation, and output bounds.
  DPE arithmetic is unchanged; this binary deliberately does not enable the
  broad UBSan checks associated with the separate DPE repair.

The production patch changes only the two upstream algorithm fragments and
their byte counts/SHA-256 records in each backend manifest. Fixed-precision
capacities, packed ML2 state, workspace allocation, and retry policy are
unchanged. The regression sources and runner were reviewed but not compiled
or executed when this port was prepared. No host-side tests, ARM
cross-compilation, or KeyGen/Sign/Verify run was performed for this change.
