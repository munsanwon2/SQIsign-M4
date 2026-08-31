/* Genuine compact-SQIsign/Qlapoti library bindings for the pqm4 adapter.
 * The curated backend compile flags select lvl1, lvl3, or lvl5, so the
 * namespace macros below resolve to the matching fixed-precision symbols. */
#include "sqisign_m4_backend.h"

#include <encoded_sizes.h>
#include <sig.h>
#include <signature.h>
#include <verification.h>

#include <string.h>

/*
 * The STM32F407 has a separately addressable 16 KiB SRAM2 bank.  Keeping the
 * decoded key/signature objects there removes their roughly 12 KiB lvl-V
 * footprint from the already tight main-stack budget.  Other targets retain
 * automatic storage, which keeps the portable/host implementation reentrant.
 */
typedef union sqisign_m4_binding_state {
    struct {
        secret_key_t secret_key;
        public_key_t public_key;
    } keypair;
    struct {
        secret_key_t secret_key;
        public_key_t public_key;
        signature_t signature;
    } sign;
    struct {
        public_key_t public_key;
        signature_t signature;
    } verify;
} sqisign_m4_binding_state_t;

#if defined(SQISIGN_M4_F407_BANKED_RAM)
_Static_assert(sizeof(sqisign_m4_binding_state_t) <= 16u * 1024u,
               "SQIsign binding state no longer fits STM32F407 SRAM2");

static sqisign_m4_binding_state_t sqisign_m4_binding_state
    __attribute__((section(".ram2.sqisign_binding"), aligned(8)));
static unsigned int sqisign_m4_binding_busy;

static sqisign_m4_binding_state_t *
sqisign_m4_binding_acquire(sqisign_m4_binding_state_t *automatic_state)
{
    (void)automatic_state;
    if (__atomic_exchange_n(&sqisign_m4_binding_busy,
                            1u,
                            __ATOMIC_ACQUIRE) != 0u)
        return NULL;
    return &sqisign_m4_binding_state;
}
#else
static sqisign_m4_binding_state_t *
sqisign_m4_binding_acquire(sqisign_m4_binding_state_t *automatic_state)
{
    return automatic_state;
}
#endif

static void
sqisign_m4_binding_release(sqisign_m4_binding_state_t *state)
{
    /* Keep decoded secret material out of persistent RAM between operations. */
    volatile unsigned char *p = (volatile unsigned char *)state;
    size_t n = sizeof(*state);
    while (n-- != 0u)
        *p++ = 0u;
#if defined(SQISIGN_M4_F407_BANKED_RAM)
    __atomic_store_n(&sqisign_m4_binding_busy, 0u, __ATOMIC_RELEASE);
#endif
}

int
sqisign_m4_backend_keypair(unsigned char *pk, unsigned char *sk)
{
    int status = -1;
#if !defined(SQISIGN_M4_F407_BANKED_RAM)
    sqisign_m4_binding_state_t automatic_state;
#endif
    sqisign_m4_binding_state_t *state;

    if (pk == NULL || sk == NULL)
        return -1;

#if defined(SQISIGN_M4_F407_BANKED_RAM)
    state = sqisign_m4_binding_acquire(NULL);
#else
    state = sqisign_m4_binding_acquire(&automatic_state);
#endif
    if (state == NULL)
        return -1;
    memset(state, 0, sizeof(*state));
    secret_key_init(&state->keypair.secret_key);

    if (!protocols_keygen(&state->keypair.public_key,
                          &state->keypair.secret_key))
        goto cleanup;
    if (!secret_key_to_bytes(sk,
                             &state->keypair.secret_key,
                             &state->keypair.public_key))
        goto cleanup;
    public_key_to_bytes(pk, &state->keypair.public_key);
    status = 0;

cleanup:
    if (status != 0) {
        memset(pk, 0, (size_t)PUBLICKEY_BYTES);
        memset(sk, 0, (size_t)SECRETKEY_BYTES);
    }
    secret_key_finalize(&state->keypair.secret_key);
    sqisign_m4_binding_release(state);
    return status;
}

int
sqisign_m4_backend_sign(unsigned char *sig,
                         size_t *siglen,
                         const unsigned char *msg,
                         size_t msglen,
                         const unsigned char *sk)
{
    int status = -1;
#if !defined(SQISIGN_M4_F407_BANKED_RAM)
    sqisign_m4_binding_state_t automatic_state;
#endif
    sqisign_m4_binding_state_t *state;

    if (sig == NULL || siglen == NULL || sk == NULL ||
        (msg == NULL && msglen != 0u)) {
        return -1;
    }
    *siglen = 0u;
#if defined(SQISIGN_M4_F407_BANKED_RAM)
    state = sqisign_m4_binding_acquire(NULL);
#else
    state = sqisign_m4_binding_acquire(&automatic_state);
#endif
    if (state == NULL)
        return -1;
    memset(state, 0, sizeof(*state));
    secret_key_init(&state->sign.secret_key);

    if (!secret_key_from_bytes(&state->sign.secret_key,
                               &state->sign.public_key,
                               sk))
        goto cleanup;
    if (!protocols_sign(&state->sign.signature,
                        &state->sign.public_key,
                        &state->sign.secret_key,
                        msg,
                        msglen))
        goto cleanup;

    signature_to_bytes(sig, &state->sign.signature);
    *siglen = (size_t)SIGNATURE_BYTES;
    status = 0;

cleanup:
    secret_key_finalize(&state->sign.secret_key);
    sqisign_m4_binding_release(state);
    return status;
}

int
sqisign_m4_backend_verify(const unsigned char *sig,
                           size_t siglen,
                           const unsigned char *msg,
                           size_t msglen,
                           const unsigned char *pk)
{
    int status = -1;
#if !defined(SQISIGN_M4_F407_BANKED_RAM)
    sqisign_m4_binding_state_t automatic_state;
#endif
    sqisign_m4_binding_state_t *state;

    if (sig == NULL || pk == NULL || (msg == NULL && msglen != 0u) ||
        siglen != (size_t)SIGNATURE_BYTES) {
        return -1;
    }

#if defined(SQISIGN_M4_F407_BANKED_RAM)
    state = sqisign_m4_binding_acquire(NULL);
#else
    state = sqisign_m4_binding_acquire(&automatic_state);
#endif
    if (state == NULL)
        return -1;
    memset(state, 0, sizeof(*state));
    if (public_key_from_bytes(&state->verify.public_key, pk) == NULL ||
        !signature_from_bytes(&state->verify.signature, sig))
        goto cleanup;

    status = protocols_verify(&state->verify.signature,
                              &state->verify.public_key,
                              msg,
                              msglen)
                 ? 0
                 : -1;

cleanup:
    sqisign_m4_binding_release(state);
    return status;
}
