# sqisign-qlapoti-lvl1 / m4

This directory is the Cortex-M4 pqm4 wrapper for Level I. It contains
the combined-message API, F407 split-RAM workspace, operation guard, entropy
bridge, and build configuration.

The wrapper is intentionally source-injection based. Run
`scripts/inject_scheme_backends.sh` to add the genuine fixed-precision
compact-SQIsign + Qlapoti export and exactly one KeyGen/Sign/Verify binding.
After injection, the audit requires the complete 9x29 field, EC, theta,
quaternion/Qlapoti, IdealToIsogeny, signature, and verification closure.

The unpopulated overlay alone remains fail-closed at link time; it is a
distribution template, not a mock implementation. The validated injector input
is the repository's `backends/lvl1` export root; its `source/` child is the
manifested payload.
