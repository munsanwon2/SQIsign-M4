
#ifndef TOOLS_H
#define TOOLS_H

#if !defined(TARGET_ARM)
#include <time.h>
#endif

// Debug printing:
// https://stackoverflow.com/questions/1644868/define-macro-for-debug-printing-in-c
#ifndef __FILE_NAME__
#define __FILE_NAME__ "NA"
#endif

#ifndef __LINE__
#define __LINE__ 0
#endif

#ifndef __func__
#define __func__ "NA"
#endif

#if defined(TARGET_ARM)
/* The production Cortex-M4 backend has no console and must never acquire a
 * libc I/O dependency through diagnostics. */
#define debug_print(message) do { (void)sizeof(message); } while (0)
#else
void sqisign_host_debug_warning(const char *message,
                                const char *file,
                                int line,
                                const char *function);

#ifndef NDEBUG
#define debug_print(message)                                                                       \
    sqisign_host_debug_warning((message), __FILE_NAME__, __LINE__, __func__)
#else
#define debug_print(message) do { (void)sizeof(message); } while (0)
#endif

clock_t tic(void);
float tac(void);                             /* time in ms since last tic */
float TAC(const char *str);                  /* same, but prints it with label 'str' */
float toc(const clock_t t);                  /* time in ms since t */
float TOC(const clock_t t, const char *str); /* same, but prints it with label 'str' */
float TOC_clock(const clock_t t, const char *str);

clock_t dclock(const clock_t t); // return the clock cycle diff between now and t
float clock_to_time(const clock_t t,
                    const char *str); // convert the number of clock cycles t to time
float clock_print(const clock_t t, const char *str);
#endif
#endif
