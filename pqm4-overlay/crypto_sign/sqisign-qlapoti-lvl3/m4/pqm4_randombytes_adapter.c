/* Bridge SQIsign's named size_t-width entropy hook to pqm4's platform RNG.
 * Include pqm4's declaration instead of redeclaring a raw symbol: the locked
 * pqm4 tree namespaces this API to PQCLEAN_randombytes via randombytes.h. */
#include "randombytes.h"

int sqisign_platform_randombytes(unsigned char *output, size_t output_length)
{
    return randombytes(output, output_length);
}
