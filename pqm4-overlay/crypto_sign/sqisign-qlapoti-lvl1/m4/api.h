#ifndef SQISIGN_QLAPOTI_LVL1_M4_API_H
#define SQISIGN_QLAPOTI_LVL1_M4_API_H

#include <stddef.h>

#define CRYPTO_ALGNAME "SQIsign-Qlapoti-L1"
#define CRYPTO_PUBLICKEYBYTES 65u
#define CRYPTO_SECRETKEYBYTES 353u
#define CRYPTO_BYTES 148u

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);
int crypto_sign(unsigned char *sm, size_t *smlen,
                const unsigned char *msg, size_t len,
                const unsigned char *sk);
int crypto_sign_open(unsigned char *m, size_t *mlen,
                     const unsigned char *sm, size_t smlen,
                     const unsigned char *pk);

#endif
