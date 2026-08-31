#ifndef SQISIGN_M4_OPERATION_GUARD_H
#define SQISIGN_M4_OPERATION_GUARD_H

#include "qlapoti_workspace.h"

typedef int (*sqisign_m4_core_fn)(void *context);

enum {
    SQISIGN_M4_GUARD_OK = 0,
    SQISIGN_M4_GUARD_BAD_ARGUMENT = -1,
    SQISIGN_M4_GUARD_CORE_FAILED = -2,
    SQISIGN_M4_GUARD_WORKSPACE_FAILED = -3,
    SQISIGN_M4_GUARD_BUSY = -4
};

extern volatile qlw_stats_t sqisign_m4_last_workspace_stats;
extern volatile int sqisign_m4_last_guard_status;

int sqisign_m4_guard_run(sqisign_m4_core_fn core, void *context);

#endif
