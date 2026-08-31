/* Offline-generated D:E=2:1 theta schedule.
 * Baseline theta_isogenies.c SHA-256:
 * 4748b38c2e39abeeda4a1cf87b19a9ce376d8aba0f61367c6c3a5d9458a3eb7b
 * Split payload SHA-256:
 * 4a281b5a4ee012bea165799c95185acf025a9d256ed3123f1aa590724b3a6aa0
 * Generator manifest SHA-256:
 * eee5dccba4b87bd75c6f8fe9f0940a8013c22a25d5637aeaf1d7f1f6084e3a84 */
#ifndef SQISIGN_THETA_SCHEDULE_DE2_H
#define SQISIGN_THETA_SCHEDULE_DE2_H

#include <stddef.h>
#include <stdint.h>

enum {
    THETA_SCHEDULE_DE2_MAX_SLOTS = 9,
    THETA_SCHEDULE_DE2_TABLE_BYTES = 674,
    THETA_SCHEDULE_DE2_CONST_BYTES = 690
};

/* Row d=3..9 begins at offsets[d-3]; the final entry is the end. */
static const uint16_t theta_schedule_de2_offsets[8] = {
    0, 9, 29, 68, 139, 252, 419, 674
};

static const uint8_t theta_schedule_de2_splits[674] = {
    1, 1, 2, 2, 3, 4, 4, 5, 6, 1, 1, 1, 2, 2, 3, 3, 3, 4, 5, 5,
    6, 7, 7, 7, 7, 8, 9, 10, 11, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5,
    5, 5, 5, 5, 6, 7, 8, 8, 8, 8, 9, 10, 11, 12, 12, 12, 12, 12, 12, 12,
    12, 13, 14, 15, 16, 17, 18, 19, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5,
    5, 5, 5, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 9, 10, 11, 12, 13, 13,
    13, 13, 13, 13, 13, 13, 14, 15, 16, 17, 18, 19, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 1,
    1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 6, 7, 8, 8,
    8, 8, 8, 8, 8, 8, 9, 10, 11, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 14, 15, 16, 17, 18, 19, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
    21, 21, 21, 21, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 33, 33, 33,
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33,
    33, 33, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 1, 1, 1, 2, 2, 2, 3, 3,
    3, 3, 4, 5, 5, 5, 5, 5, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    9, 10, 11, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
    21, 21, 21, 21, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 34, 34,
    34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
    34, 34, 34, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
    54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 1,
    1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 6, 7, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 9, 10, 11, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 14, 15, 16, 17, 18, 19, 20, 21, 21, 21, 21, 21, 21, 21, 21,
    21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22, 23, 24, 25, 26, 27, 28,
    29, 30, 31, 32, 33, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
    34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 35,
    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
    55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55,
    55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55,
    55, 55, 55, 55, 55, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
    70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 88,
    88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88,
    88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88,
};

_Static_assert(sizeof(theta_schedule_de2_offsets) == 16,
               "theta schedule offset bytes changed");
_Static_assert(sizeof(theta_schedule_de2_splits) == 674,
               "theta schedule split bytes changed");

/* Fail closed: there is deliberately no half-split fallback. */
static inline int
theta_schedule_de2_lookup(unsigned current, uint16_t remaining, uint16_t *num_dbls)
{
    if (num_dbls == NULL || remaining < 2 || current > 7)
        return 0;
    const unsigned depth = THETA_SCHEDULE_DE2_MAX_SLOTS - current;
    if (depth == 2) {
        if (remaining > 4)
            return 0;
        *num_dbls = remaining - 1;
        return 1;
    }
    if (depth < 3 || depth > THETA_SCHEDULE_DE2_MAX_SLOTS)
        return 0;
    const unsigned row = depth - 3;
    const uint16_t begin = theta_schedule_de2_offsets[row];
    const uint16_t end = theta_schedule_de2_offsets[row + 1];
    const uint16_t maximum = (uint16_t)(end - begin + 1);
    if (remaining > maximum)
        return 0;
    const uint16_t split = theta_schedule_de2_splits[begin + remaining - 2];
    if (split == 0 || split >= remaining)
        return 0;
    *num_dbls = split;
    return 1;
}

#endif

