#include "api.h"

#include <stddef.h>
#include <stdint.h>

#include "m4_operation_guard.h"
#include "secure_bzero.h"
#include "sqisign_m4_backend.h"

static int range_end(const void *ptr, size_t len, uintptr_t *begin, uintptr_t *end) {
    uintptr_t b;
    if (ptr == 0 && len != 0u) {
        return -1;
    }
    b = (uintptr_t)ptr;
    if (len > (size_t)(UINTPTR_MAX - b)) {
        return -1;
    }
    *begin = b;
    *end = b + (uintptr_t)len;
    return 0;
}

static int ranges_overlap(const void *a, size_t alen, const void *b, size_t blen) {
    uintptr_t ab, ae, bb, be;
    if (alen == 0u || blen == 0u) {
        return 0;
    }
    if (range_end(a, alen, &ab, &ae) != 0 || range_end(b, blen, &bb, &be) != 0) {
        return 1;
    }
    return ab < be && bb < ae;
}

static void copy_forward(unsigned char *dst, const unsigned char *src, size_t len) {
    size_t i;
    for (i = 0u; i < len; ++i) {
        dst[i] = src[i];
    }
}

static void copy_backward(unsigned char *dst, const unsigned char *src, size_t len) {
    while (len != 0u) {
        --len;
        dst[len] = src[len];
    }
}

typedef struct { unsigned char *pk; unsigned char *sk; } keypair_context_t;
static int keypair_core(void *opaque) {
    keypair_context_t *ctx = (keypair_context_t *)opaque;
    return sqisign_m4_backend_keypair(ctx->pk, ctx->sk);
}

typedef struct {
    unsigned char *sig;
    size_t siglen;
    const unsigned char *msg;
    size_t msglen;
    const unsigned char *sk;
} sign_context_t;
static int sign_core(void *opaque) {
    sign_context_t *ctx = (sign_context_t *)opaque;
    ctx->siglen = 0u;
    return sqisign_m4_backend_sign(ctx->sig, &ctx->siglen,
                                  ctx->msg, ctx->msglen, ctx->sk);
}

typedef struct {
    const unsigned char *sig;
    size_t siglen;
    const unsigned char *msg;
    size_t msglen;
    const unsigned char *pk;
} verify_context_t;
static int verify_core(void *opaque) {
    verify_context_t *ctx = (verify_context_t *)opaque;
    return sqisign_m4_backend_verify(ctx->sig, ctx->siglen,
                                     ctx->msg, ctx->msglen, ctx->pk);
}

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
    keypair_context_t ctx;
    int status;
    if (pk == 0 || sk == 0 ||
        ranges_overlap(pk, CRYPTO_PUBLICKEYBYTES, sk, CRYPTO_SECRETKEYBYTES)) {
        return -1;
    }
    ctx.pk = pk;
    ctx.sk = sk;
    status = sqisign_m4_guard_run(keypair_core, &ctx);
    if (status != SQISIGN_M4_GUARD_OK) {
        sqisign_m4_secure_bzero(pk, CRYPTO_PUBLICKEYBYTES);
        sqisign_m4_secure_bzero(sk, CRYPTO_SECRETKEYBYTES);
        return -1;
    }
    return 0;
}

int crypto_sign(unsigned char *sm, size_t *smlen,
                const unsigned char *msg, size_t len,
                const unsigned char *sk) {
    sign_context_t ctx;
    const unsigned char *message_for_backend = msg;
    size_t combined_len;
    int status;

    if (smlen == 0) {
        return -1;
    }
    *smlen = 0u;
    if (sm == 0 || sk == 0 || (msg == 0 && len != 0u) ||
        len > SIZE_MAX - (size_t)CRYPTO_BYTES) {
        return -1;
    }
    combined_len = (size_t)CRYPTO_BYTES + len;
    if (sm != msg && ranges_overlap(sm, combined_len, msg, len)) {
        return -1;
    }

    /* Exact in-place operation is supported, as expected by PQClean-style APIs. */
    if (sm == msg && len != 0u) {
        copy_backward(sm + CRYPTO_BYTES, sm, len);
        message_for_backend = sm + CRYPTO_BYTES;
    }

    ctx.sig = sm;
    ctx.siglen = 0u;
    ctx.msg = message_for_backend;
    ctx.msglen = len;
    ctx.sk = sk;
    status = sqisign_m4_guard_run(sign_core, &ctx);
    if (status != SQISIGN_M4_GUARD_OK || ctx.siglen != (size_t)CRYPTO_BYTES) {
        sqisign_m4_secure_bzero(sm, combined_len);
        return -1;
    }
    if (sm != msg && len != 0u) {
        copy_forward(sm + CRYPTO_BYTES, msg, len);
    }
    *smlen = combined_len;
    return 0;
}

int crypto_sign_open(unsigned char *m, size_t *mlen,
                     const unsigned char *sm, size_t smlen,
                     const unsigned char *pk) {
    size_t message_len;
    verify_context_t ctx;
    int status;

    if (mlen == 0) {
        return -1;
    }
    *mlen = 0u;
    if (m == 0 || sm == 0 || pk == 0 || smlen < (size_t)CRYPTO_BYTES) {
        return -1;
    }
    message_len = smlen - (size_t)CRYPTO_BYTES;
    if (m != sm && ranges_overlap(m, message_len, sm, smlen)) {
        return -1;
    }
    ctx.sig = sm;
    ctx.siglen = CRYPTO_BYTES;
    ctx.msg = sm + CRYPTO_BYTES;
    ctx.msglen = message_len;
    ctx.pk = pk;
    status = sqisign_m4_guard_run(verify_core, &ctx);
    if (status != SQISIGN_M4_GUARD_OK) {
        sqisign_m4_secure_bzero(m, message_len);
        return -1;
    }
    if (message_len != 0u) {
        copy_forward(m, sm + CRYPTO_BYTES, message_len);
    }
    *mlen = message_len;
    return 0;
}
