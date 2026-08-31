#include "m4_operation_guard.h"

volatile qlw_stats_t sqisign_m4_last_workspace_stats;
volatile int sqisign_m4_last_guard_status;
static volatile unsigned sqisign_m4_guard_active;

int sqisign_m4_guard_run(sqisign_m4_core_fn core, void *context) {
    int core_status;
    int result;
    unsigned workspace_status;
    qlw_stats_t stats;

    if (core == 0) {
        sqisign_m4_last_guard_status = SQISIGN_M4_GUARD_BAD_ARGUMENT;
        return SQISIGN_M4_GUARD_BAD_ARGUMENT;
    }
    if (__atomic_exchange_n(&sqisign_m4_guard_active, 1u,
                            __ATOMIC_ACQUIRE) != 0u) {
        sqisign_m4_last_guard_status = SQISIGN_M4_GUARD_BUSY;
        return SQISIGN_M4_GUARD_BUSY;
    }

    qlw_operation_begin();
    core_status = core(context);
    workspace_status = qlw_operation_end();
    stats = qlw_get_stats();
    sqisign_m4_last_workspace_stats = stats;
    qlw_reset();

    if (core_status != 0) {
        result = SQISIGN_M4_GUARD_CORE_FAILED;
    } else if (workspace_status != QLW_STATUS_OK) {
        result = SQISIGN_M4_GUARD_WORKSPACE_FAILED;
    } else {
        result = SQISIGN_M4_GUARD_OK;
    }
    sqisign_m4_last_guard_status = result;
    __atomic_store_n(&sqisign_m4_guard_active, 0u, __ATOMIC_RELEASE);
    return result;
}
