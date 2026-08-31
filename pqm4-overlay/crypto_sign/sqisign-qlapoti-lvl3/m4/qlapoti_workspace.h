#ifndef QLAPOTI_WORKSPACE_H
#define QLAPOTI_WORKSPACE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef QLW_WORKSPACE_BYTES
/* 25,120-byte persistent signing cap + 30,528-byte largest phase cap,
 * allocator headers, and 5 KiB of fail-closed headroom. */
#define QLW_WORKSPACE_BYTES ((size_t)61440u)
#endif

#ifndef QLW_SECTION_NAME
#define QLW_SECTION_NAME ".bss.qlapoti_workspace"
#endif

enum {
    QLW_STATUS_OK = 0u,
    QLW_STATUS_OOM = 1u << 0,
    QLW_STATUS_CORRUPT = 1u << 1,
    QLW_STATUS_LEAK = 1u << 2
};

typedef struct qlw_stats {
    size_t capacity_bytes;
    size_t current_bytes;
    size_t peak_bytes;
    size_t live_blocks;
    unsigned status_flags;
} qlw_stats_t;

/* Start a fresh signing/key-generation operation. This securely clears the
 * whole fixed workspace and resets all error flags and peak counters. */
void qlw_operation_begin(void);

/* Return QLW_STATUS_* flags. A nonzero result must make the surrounding
 * cryptographic operation fail. This function does not erase the arena;
 * call qlw_reset() after copying any desired diagnostics. */
unsigned qlw_operation_end(void);

/* Securely erase and reinitialize the arena. */
void qlw_reset(void);

void *qlw_malloc(size_t size);
void *qlw_calloc(size_t count, size_t size);
void *qlw_realloc(void *ptr, size_t size);
void qlw_free(void *ptr);

/* Fixed-backend scratch hooks consumed by sqisign_m4_workspace.h. */
void *sqisign_m4_workspace_acquire(size_t bytes);
void sqisign_m4_workspace_release(void *ptr);

qlw_stats_t qlw_get_stats(void);

/* Exposed for linker-map and board diagnostics. Do not write through these
 * pointers outside this allocator. */
extern unsigned char qlw_workspace_storage[QLW_WORKSPACE_BYTES];
const unsigned char *qlw_workspace_begin(void);
const unsigned char *qlw_workspace_end(void);

#ifdef __cplusplus
}
#endif

#endif
