# SQIsign-Qlapoti Level I

The target Cortex-M4 implementation is in `m4/`. The repository installer
injects the content-addressed Level-I export from `backends/lvl1`. That export
root holds `library-manifest.json`; its `source/` child is the manifested
payload.

The portable companion is installed under
`mupq/crypto_sign/sqisign-qlapoti-lvl1/ref/`. Target and reference use the
same algorithm and wire format.
