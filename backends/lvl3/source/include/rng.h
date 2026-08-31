// SPDX-License-Identifier: Apache-2.0

#ifndef rng_h
#define rng_h

#include <stddef.h>
#include <sqisign_namespace.h>

#if defined(TARGET_ARM) || defined(SQISIGN_USE_PLATFORM_RANDOMBYTES)
/**
 * Platform entropy provider supplied by the board or harness integration.
 *
 * The implementation must fill exactly output_length bytes with
 * cryptographically secure random data and return 0 on success. Any nonzero
 * return value reports failure and must leave the caller treating the output
 * as unusable. No weak/default implementation is provided. Bare-metal builds
 * select this contract through TARGET_ARM. A host harness that deliberately
 * supplies a platform-style provider may opt in narrowly with
 * SQISIGN_USE_PLATFORM_RANDOMBYTES, without enabling any other ARM ABI or
 * memory-layout behavior.
 */
int sqisign_platform_randombytes(unsigned char *output, size_t output_length);
#endif

/**
 * Randombytes initialization.
 * Initialization may be needed for some random number generators (e.g. CTR-DRBG).
 *
 * @param[in] entropy_input 48 bytes entropy input
 * @param[in] personalization_string Personalization string
 * @param[in] security_strength Security string
 */
SQISIGN_API
void randombytes_init(unsigned char *entropy_input,
                      unsigned char *personalization_string,
                      int security_strength);

/**
 * Random byte generation using /dev/urandom.
 * The caller is responsible to allocate sufficient memory to hold x.
 *
 * @param[out] x Memory to hold the random bytes.
 * @param[in] xlen Number of random bytes to be generated
 * @return int 0 on success, -1 otherwise
 */
SQISIGN_API
int randombytes_select(unsigned char *x, unsigned long long xlen);

/**
 * Random byte generation.
 * The caller is responsible to allocate sufficient memory to hold x.
 *
 * @param[out] x Memory to hold the random bytes.
 * @param[in] xlen Number of random bytes to be generated
 * @return int 0 on success, -1 otherwise
 */
SQISIGN_API
int randombytes(unsigned char *x, unsigned long long xlen);

/**
 * Production-library random byte request with a target-native length ABI.
 *
 * pqm4 exposes ``randombytes(uint8_t *, size_t)``.  Calling that symbol
 * through SQIsign's historical ``unsigned long long`` prototype is not ABI
 * compatible on 32-bit ARM (AAPCS passes the length in different registers).
 * Keep the public host API unchanged, but make embedded production call sites
 * and explicitly opted-in harnesses use this size_t-width wrapper and the
 * named platform hook.
 */
static inline int
sqisign_randombytes(unsigned char *x, size_t xlen)
{
#if defined(TARGET_ARM) || defined(SQISIGN_USE_PLATFORM_RANDOMBYTES)
    return sqisign_platform_randombytes(x, xlen);
#else
    return randombytes(x, (unsigned long long)xlen);
#endif
}

/* Wipe temporary entropy buffers through volatile stores so the compiler
 * cannot elide the cleanup after their contents have been consumed. */
static inline void
sqisign_randombytes_wipe(void *buffer, size_t length)
{
    volatile unsigned char *cursor = (volatile unsigned char *)buffer;
    while (length-- != 0u)
        *cursor++ = 0u;
}

#endif /* rng_h */
