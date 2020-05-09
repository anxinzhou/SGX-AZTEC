//
// Created by anxin on 2020/5/6.
//

#ifndef PROOF_UTIL_H
#define PROOF_UTIL_H
#include <stdint.h>


extern unsigned long next;
extern const unsigned int MY_RAND_MAX;

/* RAND_MAX assumed to be 32767 */
static inline int myrand_int() {
    next = next * 1103515245 + 12345;
    return((unsigned)(next/65536) % 32768);
}

static inline double myrandom() {
    return (double) myrand_int()/MY_RAND_MAX;
}

static inline void mysrand(unsigned seed) {
    next = seed;
}

static inline uint32_t oblivious_assign_CMOV(uint8_t pred, uint32_t t_val, uint32_t f_val) {
    return pred ? t_val:f_val;
}
#endif //PROOF_UTIL_H
