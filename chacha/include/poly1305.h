#pragma once

#include "array_interface.h"

#include <iostream>
#include <gmpxx.h>

#define YKM_POLY1305_SIZE 3

/*
p = fffffffffffffffb ffffffffffffffff 3
*/

void poly1305_mac(uint8_t*, uint8_t*, size_t, uint8_t*);
void poly1305_key_gen(uint8_t*, uint8_t*, uint8_t*);

inline void DivMod(mp_limb_t *r, mp_limb_t *x, mp_limb_t *y) {
    mp_limb_t p[YKM_POLY1305_SIZE] = {0xfffffffffffffffb, 0xffffffffffffffff, 0x3};
    mp_limb_t xy[2*YKM_POLY1305_SIZE];
    mp_limb_t q[YKM_POLY1305_SIZE + 1];

    mpn_mul_n(xy, x, y, YKM_POLY1305_SIZE);
#ifdef NDEBUG
    printf("%lx ", xy[3]);
    printf("%lx ", xy[2]);
    printf("%lx ", xy[1]);
    printf("%lx\n",xy[0]);
    puts("");
#endif

    mpn_tdiv_qr(q, r, 0, xy, 2*YKM_POLY1305_SIZE, p, YKM_POLY1305_SIZE);
}

inline void add_3(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y) {
    mpn_add_n(z, (const mp_limb_t*)x, (const mp_limb_t*)y, YKM_POLY1305_SIZE);
}

inline void poly1305aes_test_clamp(uint8_t r[16]) {
    r[3] &= 15;
    r[7] &= 15;
    r[11] &= 15;
    r[15] &= 15;
    r[4] &= 252;
    r[8] &= 252;
    r[12] &= 252;
}
