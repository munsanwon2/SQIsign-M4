#ifndef SQISIGN_M4_WORKSPACE_H
#define SQISIGN_M4_WORKSPACE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Bounded transient storage supplied by the Cortex-M4 integration layer.
 * The pqm4 port maps this allocator onto its single CCM NOLOAD arena.  These
 * hooks deliberately are not libc malloc/free and are unavailable to host
 * builds, which keep ordinary automatic storage for concurrency tests. */
void *sqisign_m4_workspace_acquire(size_t bytes);
void sqisign_m4_workspace_release(void *ptr);

#ifdef __cplusplus
}
#endif

#endif
