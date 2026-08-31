/*
 * On ARM, route hosted-style transient allocation used by the imported
 * backend into the bounded workspace.  Host pqm4 testvector executables keep
 * their native allocator; the real backend should nevertheless be tested with
 * the arena in its own native test configuration.
 */
#include "qlapoti_workspace.h"
#include <stddef.h>

#if defined(__arm__) || defined(__thumb__) || defined(QLAPOTI_OVERRIDE_LIBC_ALLOC)
void *malloc(size_t size) { return qlw_malloc(size); }
void *calloc(size_t count, size_t size) { return qlw_calloc(count, size); }
void *realloc(void *ptr, size_t size) { return qlw_realloc(ptr, size); }
void free(void *ptr) { qlw_free(ptr); }
#endif
