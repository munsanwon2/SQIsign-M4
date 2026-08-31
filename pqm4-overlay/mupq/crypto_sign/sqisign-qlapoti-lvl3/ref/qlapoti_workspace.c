#include "qlapoti_workspace.h"

#include <stdalign.h>
#include <stdbool.h>
#include <stdint.h>

#define QLW_MAGIC UINT32_C(0x514c5742) /* "QLWB" */
#define QLW_USED UINT32_C(0x00000001)
#define QLW_MIN_SPLIT ((size_t)16u)

#if defined(__GNUC__) || defined(__clang__)
#define QLW_SECTION __attribute__((section(QLW_SECTION_NAME), aligned(16), used))
#else
#define QLW_SECTION
#endif

typedef struct qlw_block {
    uint32_t magic;
    uint32_t flags;
    size_t size;
    struct qlw_block *prev;
    struct qlw_block *next;
} qlw_block_t;

unsigned char qlw_workspace_storage[QLW_WORKSPACE_BYTES] QLW_SECTION;
static qlw_block_t *g_head;
static qlw_stats_t g_stats;
static bool g_initialized;

static size_t qlw_alignment(void) {
    const size_t a = alignof(max_align_t);
    return a < (size_t)8u ? (size_t)8u : a;
}

static bool add_overflow_size(size_t a, size_t b, size_t *out) {
    if (a > SIZE_MAX - b) {
        return true;
    }
    *out = a + b;
    return false;
}

static bool mul_overflow_size(size_t a, size_t b, size_t *out) {
    if (a != 0u && b > SIZE_MAX / a) {
        return true;
    }
    *out = a * b;
    return false;
}

static bool align_up_size(size_t n, size_t *out) {
    const size_t a = qlw_alignment();
    const size_t mask = a - 1u;
    size_t tmp;
    if ((a & mask) != 0u || add_overflow_size(n, mask, &tmp)) {
        return false;
    }
    *out = tmp & ~mask;
    return true;
}

static size_t header_bytes(void) {
    size_t out = 0u;
    (void)align_up_size(sizeof(qlw_block_t), &out);
    return out;
}

static unsigned char *block_payload(qlw_block_t *block) {
    return ((unsigned char *)block) + header_bytes();
}

static const unsigned char *block_payload_const(const qlw_block_t *block) {
    return ((const unsigned char *)block) + header_bytes();
}

static void secure_zero(void *ptr, size_t len) {
    volatile unsigned char *p = (volatile unsigned char *)ptr;
    while (len-- != 0u) {
        *p++ = 0u;
    }
}

static void copy_bytes(void *dst, const void *src, size_t len) {
    unsigned char *d = (unsigned char *)dst;
    const unsigned char *s = (const unsigned char *)src;
    while (len-- != 0u) {
        *d++ = *s++;
    }
}

static bool span_in_workspace(const void *ptr, size_t len) {
    const uintptr_t lo = (uintptr_t)&qlw_workspace_storage[0];
    const uintptr_t hi = lo + (uintptr_t)sizeof(qlw_workspace_storage);
    const uintptr_t p = (uintptr_t)ptr;
    if (p < lo || p > hi) {
        return false;
    }
    return len <= (size_t)(hi - p);
}

static void set_corrupt(void) {
    g_stats.status_flags |= QLW_STATUS_CORRUPT;
}

static void initialize_if_needed(void) {
    size_t hdr;
    if (g_initialized) {
        return;
    }
    secure_zero(qlw_workspace_storage, sizeof(qlw_workspace_storage));
    hdr = header_bytes();
    if (hdr == 0u || hdr >= sizeof(qlw_workspace_storage)) {
        g_head = NULL;
        g_stats = (qlw_stats_t){0};
        g_stats.status_flags = QLW_STATUS_CORRUPT;
        g_initialized = true;
        return;
    }
    g_head = (qlw_block_t *)(void *)qlw_workspace_storage;
    g_head->magic = QLW_MAGIC;
    g_head->flags = 0u;
    g_head->size = sizeof(qlw_workspace_storage) - hdr;
    g_head->prev = NULL;
    g_head->next = NULL;
    g_stats.capacity_bytes = g_head->size;
    g_stats.current_bytes = 0u;
    g_stats.peak_bytes = 0u;
    g_stats.live_blocks = 0u;
    g_stats.status_flags = QLW_STATUS_OK;
    g_initialized = true;
}

static bool validate_chain(void) {
    const size_t hdr = header_bytes();
    const size_t max_blocks = sizeof(qlw_workspace_storage) / (hdr == 0u ? 1u : hdr) + 1u;
    const unsigned char *arena_end = qlw_workspace_storage + sizeof(qlw_workspace_storage);
    const qlw_block_t *prev = NULL;
    const qlw_block_t *cur;
    size_t count = 0u;

    initialize_if_needed();
    if (g_head == NULL) {
        set_corrupt();
        return false;
    }
    cur = g_head;
    while (cur != NULL) {
        const unsigned char *payload;
        const unsigned char *end;
        if (++count > max_blocks || !span_in_workspace(cur, hdr)) {
            set_corrupt();
            return false;
        }
        if (cur->magic != QLW_MAGIC || cur->prev != prev ||
            (cur->flags & ~QLW_USED) != 0u) {
            set_corrupt();
            return false;
        }
        payload = block_payload_const(cur);
        if (!span_in_workspace(payload, cur->size)) {
            set_corrupt();
            return false;
        }
        end = payload + cur->size;
        if (cur->next != NULL) {
            if ((const unsigned char *)cur->next != end || end >= arena_end) {
                set_corrupt();
                return false;
            }
        } else if (end != arena_end) {
            set_corrupt();
            return false;
        }
        prev = cur;
        cur = cur->next;
    }
    return true;
}

static qlw_block_t *find_block(void *ptr) {
    qlw_block_t *cur;
    if (ptr == NULL || !validate_chain()) {
        return NULL;
    }
    for (cur = g_head; cur != NULL; cur = cur->next) {
        if ((void *)block_payload(cur) == ptr) {
            return cur;
        }
    }
    return NULL;
}

static void split_block(qlw_block_t *block, size_t wanted) {
    const size_t hdr = header_bytes();
    const size_t remaining = block->size - wanted;
    qlw_block_t *tail;
    if (remaining < hdr + QLW_MIN_SPLIT) {
        return;
    }
    tail = (qlw_block_t *)(void *)(block_payload(block) + wanted);
    tail->magic = QLW_MAGIC;
    tail->flags = 0u;
    tail->size = remaining - hdr;
    tail->prev = block;
    tail->next = block->next;
    if (tail->next != NULL) {
        tail->next->prev = tail;
    }
    block->next = tail;
    block->size = wanted;
}

static void merge_with_next(qlw_block_t *block) {
    qlw_block_t *next = block->next;
    if (next == NULL || (next->flags & QLW_USED) != 0u) {
        return;
    }
    block->size += header_bytes() + next->size;
    block->next = next->next;
    if (block->next != NULL) {
        block->next->prev = block;
    }
    secure_zero(next, header_bytes());
}

void qlw_operation_begin(void) {
    qlw_reset();
}

unsigned qlw_operation_end(void) {
    initialize_if_needed();
    if (!validate_chain()) {
        g_stats.status_flags |= QLW_STATUS_CORRUPT;
    }
    if (g_stats.live_blocks != 0u || g_stats.current_bytes != 0u) {
        g_stats.status_flags |= QLW_STATUS_LEAK;
    }
    return g_stats.status_flags;
}

void qlw_reset(void) {
    secure_zero(qlw_workspace_storage, sizeof(qlw_workspace_storage));
    g_head = NULL;
    g_stats = (qlw_stats_t){0};
    g_initialized = false;
    initialize_if_needed();
}

void *qlw_malloc(size_t size) {
    qlw_block_t *cur;
    size_t wanted;
    initialize_if_needed();
    if (size == 0u) {
        size = 1u;
    }
    if (!align_up_size(size, &wanted) || !validate_chain()) {
        g_stats.status_flags |= QLW_STATUS_OOM;
        return NULL;
    }
    for (cur = g_head; cur != NULL; cur = cur->next) {
        if ((cur->flags & QLW_USED) == 0u && cur->size >= wanted) {
            split_block(cur, wanted);
            cur->flags |= QLW_USED;
            g_stats.current_bytes += cur->size;
            g_stats.live_blocks += 1u;
            if (g_stats.current_bytes > g_stats.peak_bytes) {
                g_stats.peak_bytes = g_stats.current_bytes;
            }
            return block_payload(cur);
        }
    }
    g_stats.status_flags |= QLW_STATUS_OOM;
    return NULL;
}

void *qlw_calloc(size_t count, size_t size) {
    size_t total;
    void *ptr;
    if (mul_overflow_size(count, size, &total)) {
        initialize_if_needed();
        g_stats.status_flags |= QLW_STATUS_OOM;
        return NULL;
    }
    ptr = qlw_malloc(total);
    if (ptr != NULL) {
        secure_zero(ptr, total);
    }
    return ptr;
}

void qlw_free(void *ptr) {
    qlw_block_t *block;
    if (ptr == NULL) {
        return;
    }
    block = find_block(ptr);
    if (block == NULL || (block->flags & QLW_USED) == 0u) {
        set_corrupt();
        return;
    }
    if (g_stats.current_bytes < block->size || g_stats.live_blocks == 0u) {
        set_corrupt();
        return;
    }
    g_stats.current_bytes -= block->size;
    g_stats.live_blocks -= 1u;
    secure_zero(block_payload(block), block->size);
    block->flags &= ~QLW_USED;
    merge_with_next(block);
    if (block->prev != NULL && (block->prev->flags & QLW_USED) == 0u) {
        block = block->prev;
        merge_with_next(block);
    }
}

void *sqisign_m4_workspace_acquire(size_t bytes) {
    return qlw_malloc(bytes);
}

void sqisign_m4_workspace_release(void *ptr) {
    qlw_free(ptr);
}

void *qlw_realloc(void *ptr, size_t size) {
    qlw_block_t *block;
    size_t wanted;
    size_t old_size;
    if (ptr == NULL) {
        return qlw_malloc(size);
    }
    if (size == 0u) {
        qlw_free(ptr);
        return NULL;
    }
    if (!align_up_size(size, &wanted)) {
        initialize_if_needed();
        g_stats.status_flags |= QLW_STATUS_OOM;
        return NULL;
    }
    block = find_block(ptr);
    if (block == NULL || (block->flags & QLW_USED) == 0u) {
        set_corrupt();
        return NULL;
    }
    old_size = block->size;
    if (wanted <= old_size) {
        split_block(block, wanted);
        g_stats.current_bytes -= old_size - block->size;
        if (block->next != NULL && (block->next->flags & QLW_USED) == 0u) {
            merge_with_next(block->next);
        }
        return ptr;
    }
    if (block->next != NULL && (block->next->flags & QLW_USED) == 0u &&
        old_size + header_bytes() + block->next->size >= wanted) {
        merge_with_next(block);
        split_block(block, wanted);
        g_stats.current_bytes += block->size - old_size;
        if (g_stats.current_bytes > g_stats.peak_bytes) {
            g_stats.peak_bytes = g_stats.current_bytes;
        }
        return ptr;
    }
    {
        void *replacement = qlw_malloc(size);
        if (replacement == NULL) {
            return NULL;
        }
        copy_bytes(replacement, ptr, old_size < size ? old_size : size);
        qlw_free(ptr);
        return replacement;
    }
}

qlw_stats_t qlw_get_stats(void) {
    initialize_if_needed();
    return g_stats;
}

const unsigned char *qlw_workspace_begin(void) {
    return qlw_workspace_storage;
}

const unsigned char *qlw_workspace_end(void) {
    return qlw_workspace_storage + sizeof(qlw_workspace_storage);
}
