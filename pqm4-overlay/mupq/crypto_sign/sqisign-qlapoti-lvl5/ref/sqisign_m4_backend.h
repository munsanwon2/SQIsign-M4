#ifndef SQISIGN_M4_BACKEND_H
#define SQISIGN_M4_BACKEND_H

#include <stddef.h>
#include <stdint.h>

#include "api.h"

/*
 * Genuine-backend contract.  This archive intentionally supplies declarations
 * only.  A final link without the real compact-SQIsign + Qlapoti implementation
 * must fail; do not add fixed keys, fixed signatures, or verification-only
 * compatibility functions here.
 */
int sqisign_m4_backend_keypair(
    uint8_t *pk, uint8_t *sk);

int sqisign_m4_backend_sign(
    uint8_t *sig, size_t *siglen,
    const uint8_t *msg, size_t msglen,
    const uint8_t *sk);

int sqisign_m4_backend_verify(
    const uint8_t *sig, size_t siglen,
    const uint8_t *msg, size_t msglen,
    const uint8_t *pk);

#endif
