# Compact SQIsign parameters on Cortex-M4

The port uses the Round-2 protocol with compact-SQIsign Section 4
parameters. The host and Cortex-M4 builds use the same radix-32 field
representation and wire format.

| Level | Prime | λ | Torsion exponent f | Response exponent | Mixing degree bits | Integer limbs |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| I | \(3\cdot 2^{324}-1\) | 128 | 324 | 163 | 582 | 34 |
| III | \(27\cdot 2^{500}-1\) | 192 | 500 | 253 | 889 | 53 |
| V | \(17\cdot 2^{664}-1\) | 256 | 664 | 335 | 1181 | 69 |

Key generation and commitment sampling use the first prime strictly above
\(p 2^{2\lambda}\). The SHAKE challenge uses one hash iteration and λ
effective bits. The challenge ideal uses the full exponent f.

| Level | Field storage | Scalar words | Field wire bytes | Secret ideal integer bytes | Exponent bytes | Matrix entry bytes | PK / SK / signature bytes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| I | 12 × 28 bits | 12 | 48 | 48 | 1 | 21 | 97 / 501 / 200 |
| III | 18 × 29 bits | 16 | 64 | 66 | 1 | 32 | 129 / 711 / 284 |
| V | 24 × 28 bits | 22 | 88 | 88 | 2 | 43 | 177 / 953 / 386 |

The integer arrays contain 64-bit two's-complement limbs. The KeyGen
magnitude budgets are 1965/3040/4020 bits and the signing budgets are
2169/3334/4395 bits. Allocated signed magnitude capacities after limb rounding
are 2175/3391/4415 bits. Lattice operations check their intermediate width
before using the fixed-width representation.

Qlapoti uses the precomputed special order \(O_0\), its curve and action,
and its unit connecting ideal. Its initial theta-chain exponents are
322/498/662, with 321/497/661 on the transformed path. The generator
coefficients are sampled from the inclusive interval \([-2^{17},2^{17}]\).
The 20,000-iteration setting is the number of outer proposals per block.
The existing 20-block multiplier and 4,096-attempt inner caps bound the
full search. Exhausting a search produces a retry status. An arithmetic
or allocation failure produces a fatal status.
The split-prime list is generated from \(q=1\), using 2 and the first 100
primes congruent to 1 modulo 4. The level-independent list has 101 entries.

## Fixed-precision intermediates

For ML2, let \(B\) be the largest squared input numerator norm in
`diag(1,1,p,p)`. The admission check reserves nine magnitude bits above
\(B\) for \(512B\). It uses the exact diagonal near the capacity limit.
On-demand Gram entries use checked weighted products and accumulations.
An exact row update computes all four new coordinates before publishing
them and recomputing Gram entries. The ordered product path checks the
four terms of each quaternion coordinate and the common denominator.
The same-order intersection checks the trace pairings, norm product and
scaled generators before handing the eight rows to ML2.

For a Qlapoti generator proposal with reduced numerator basis columns
\(b_i\), write \(B_s=\max_i\lVert b_i\rVert_p^2\). The weighted triangle
inequality bounds the sampled numerator by \(2^{19}\sqrt{B_s}\), so its
squared norm is at most \(2^{38}B_s\) before denominator cancellation.
Within one search block the successive positive multipliers are at most
\(1+20000\cdot4096<2^{27}\). The candidate expression
\(2\lambda^2\operatorname{nrd}(\alpha_0)\) therefore needs at most 55
additional magnitude bits above the generator norm. The CVP proceeds
only after exact division by the small ideal norm and a positive result
for \(m=2^e-2\lambda^2\operatorname{nrd}(\alpha_0)/N\).
Lattice membership uses a checked adjugate and cancels denominator factors
before its final divisibility test. The checks and checked builds cover
the implementation paths in addition to the upstream Section 4 budget.

## Recreate the source exports

Use the compact-SQIsign revision in `SOURCE_URLS.lock` and the same GNU
CMake profile used to produce the manifests. The input Git clone may contain
other changes. The builder reads the locked commit through `git archive`.

```sh
python3 full-backend-tools/rebuild_compact_exports.py \
    --upstream /absolute/path/to/compact-SQIsign \
    --output /absolute/path/to/new-exports
diff -ru backends /absolute/path/to/new-exports
```

The builder verifies upstream generated-file hashes, applies the Cortex-M4
patch, selects the first entry of each generated order and curve table and
regenerates the Qlapoti constants. It then asks CMake for the complete
`sqisign_lvl1`, `sqisign_lvl3` and `sqisign_lvl5` library closures. Each
export includes the upstream parameter record beside the library manifest.
That record includes the exact prime and mixing degree, Sage seed, dependency
revisions and upstream hashes. The library manifest records the final Cortex-M4
files and hashes, including the shared field and Cornacchia include fragments.

The formulas and Sage generators are in `scripts/precomp/parameters.py` and
`scripts/precomp/precompute_sizes.sage` in the locked upstream tree. The
checked-in upstream precomputations are used by the export build.
