#ifndef SQISIGN_QLAPOTI_LVL3_M4_API_H
#define SQISIGN_QLAPOTI_LVL3_M4_API_H

#include <stddef.h>

#define CRYPTO_ALGNAME "SQIsign-Qlapoti-L3"
#define CRYPTO_PUBLICKEYBYTES 97u
#define CRYPTO_SECRETKEYBYTES 529u
#define CRYPTO_BYTES 224u

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);
int crypto_sign(unsigned char *sm, size_t *smlen,
                const unsigned char *msg, size_t len,
                const unsigned char *sk);
int crypto_sign_open(unsigned char *m, size_t *mlen,
                     const unsigned char *sm, size_t smlen,
                     const unsigned char *pk);

#endif
