#include "secure_bzero.h"

void sqisign_m4_secure_bzero(void *ptr, size_t len) {
    volatile unsigned char *p = (volatile unsigned char *)ptr;
    while (len-- != 0u) {
        *p++ = 0u;
    }
}
