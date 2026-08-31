# SQIsign-Qlapoti Level V

The target Cortex-M4 implementation is in `m4/`. The repository installer
injects the content-addressed Level-V export from `backends/lvl5`. That export
root holds `library-manifest.json`; its `source/` child is the manifested
payload.

The portable companion is installed under
`mupq/crypto_sign/sqisign-qlapoti-lvl5/ref/`. Target and reference use the
same algorithm and wire format.
