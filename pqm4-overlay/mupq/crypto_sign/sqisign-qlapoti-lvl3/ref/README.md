# sqisign-qlapoti-lvl3 / ref

This directory is the portable pqm4 reference wrapper for Level III. It
uses the same fixed-precision compact-SQIsign + Qlapoti source and 14x28
field representation as the Cortex-M4 implementation.

Run `scripts/inject_scheme_backends.sh` to add the content-addressed source
export and exactly one KeyGen/Sign/Verify binding. The reference build supplies
deterministic pqm4 testvector randomness through the same named platform hook;
it does not import the hosted OS RNG.

The unpopulated overlay alone remains fail-closed at link time. It is a
distribution template, not a mock or a claim that source injection already
occurred in this directory.
