/* Bridge SQIsign's named size_t-width entropy hook to pqm4's platform RNG. */
#if defined(TARGET_ARM) || defined(SQISIGN_USE_PLATFORM_RANDOMBYTES)
#include "randombytes.h"

int sqisign_platform_randombytes(unsigned char *output, size_t output_length)
{
    return randombytes(output, output_length);
}
#endif
