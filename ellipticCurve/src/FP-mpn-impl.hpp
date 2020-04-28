#pragma once
#include "FP.h"
#include <gmpxx.h>
#include <iostream>

mp_limb_t Fp::modulus[SIZE];

void Fp::mulInt(Fp& z, const Fp& x, int scalar) {
    mp_limb_t tmp_z[SIZE + 1] = {0};
    mp_limb_t q[2];
    mpn_copyi((mp_limb_t *)tmp_z, (const mp_limb_t *)x.value, SIZE);

    mpn_mul_1(tmp_z, (const mp_limb_t *)tmp_z, SIZE + 1, (mp_limb_t)scalar);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, SIZE + 1, (const mp_limb_t *)Fp::modulus, SIZE);
}

void Fp::setModulo(const mpz_class& v) {
    getArray(modulus, SIZE, v, v.get_mpz_t()->_mp_size);
}

void Fp::setModulo(const mp_limb_t p[SIZE]) {
    mpn_copyi((mp_limb_t *)modulus, (const mp_limb_t *)p, SIZE);
}


void Fp::neg(Fp& r, const Fp& x) {
    mpn_sub_n((mp_limb_t *)r.value, (const mp_limb_t *)Fp::modulus, (const mp_limb_t *)x.value, SIZE);
}

bool isEq(const Fp& x, const Fp& y) {
    return mpn_cmp((const mp_limb_t*)x.value, (const mp_limb_t*)y.value, SIZE) == 0;
}

void add(Fp& z, const Fp& x, const Fp& y) {
    if (mpn_add_n((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE)) {
        mp_limb_t r[SIZE];
        mpn_sub_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
        mpn_copyi((mp_limb_t *)z.value, (const mp_limb_t *)r, SIZE);
        return;
    }

    if (mpn_cmp((const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE) >= 0) {
        mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
    }
}

void add(Fp& z, const Fp& x, uint64_t scalar) {
    if (mpn_add_1((mp_limb_t *)z.value, (const mp_limb_t *)x.value, SIZE, (mp_limb_t)scalar)) {
        mp_limb_t r[SIZE];
        mpn_sub_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
        mpn_copyi((mp_limb_t *)z.value, (const mp_limb_t *)r, SIZE);
        return;
    }
    if (mpn_cmp((const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE) >= 0) {
        mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
    }
}

void sub(Fp& z, const Fp& x, const Fp& y) {
    if (mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE)) {
        mp_limb_t r[SIZE];
        mpn_add_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
        mpn_copyi((mp_limb_t *)z.value, r, SIZE);
    }
}

void mul(Fp& z, const Fp& x, const Fp& y) {
#ifdef SECP521
    mp_limb_t tmp_z[SIZE * 2];
    mp_limb_t t[SIZE*2];
    mp_limb_t s[SIZE*2];

    mpn_mul_n(tmp_z, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE);
    mod((mp_limb_t*)z.value, (const mp_limb_t *)tmp_z, (const mp_limb_t *)Fp::modulus, t, s);

#elif defined(USE_MPN)

    mp_limb_t tmp_z[SIZE * 2];
    mp_limb_t q[SIZE + 1];

    mpn_mul_n(tmp_z, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, SIZE*2, (const mp_limb_t *)Fp::modulus, SIZE);
#else
    z.value = (x.value * y.value) % Fp::modulus;
#endif
}


void invmod(Fp& r, const Fp& x) {
    mp_limb_t g[SIZE];
    mp_limb_t u[SIZE];
    mp_limb_t v[SIZE];
    mp_limb_t s[SIZE+1];
    mp_size_t sn;

    mpn_copyi(u, (const mp_limb_t*)x.value, SIZE);
    mpn_copyi(v, (const mp_limb_t*)Fp::modulus, SIZE);
    mpn_gcdext(g, s, &sn, u, SIZE, v, SIZE);
    mpn_copyi((mp_limb_t*)r.value, (const mp_limb_t*)s, SIZE);

    if (sn < 0) {
        mpn_sub_n((mp_limb_t*)r.value, (const mp_limb_t*)Fp::modulus, 
                (const mp_limb_t*)r.value, SIZE);
    }
}

void sqr(Fp &r, const Fp &x) { // r <- x^2
#ifdef SECP521
    mp_limb_t tmp_r[SIZE * 2];
    mp_limb_t t[SIZE*2];
    mp_limb_t s[SIZE*2];

    mpn_sqr(tmp_r, (const mp_limb_t *)x.value, SIZE);
    mod((mp_limb_t*)r.value, (const mp_limb_t *)tmp_r, (const mp_limb_t *)Fp::modulus, t, s);
#elif defined(USE_MPN)
    mp_limb_t tmp_r[SIZE * 2];
    mp_limb_t q[SIZE + 1];

    mpn_sqr(tmp_r, (const mp_limb_t *)x.value, SIZE);
    mpn_tdiv_qr(q, (mp_limb_t *)r.value, 0,
            tmp_r, SIZE*2, (const mp_limb_t *)Fp::modulus, SIZE);
#else
    mpz_powm_ui(r.value.get_mpz_t(), x.value.get_mpz_t(), 2, Fp::modulus.get_mpz_t());
#endif
}

bool Fp::squareRoot(Fp& r, const Fp& x) {
    mp_limb_t q[SIZE];
    mp_limb_t c[SIZE], t[SIZE], b[SIZE], z[SIZE];
    mp_limb_t tp[mpn_sec_powm_itch(SIZE, SIZE*GMP_NUMB_BITS, SIZE)];
    mp_bitcnt_t bcnt;

    mpn_copyi(t, (const mp_limb_t*)Fp::modulus, SIZE); // t = p
    mpn_rshift(t, (const mp_limb_t*)t, SIZE, 1); // t = (p-1)/2
    powMod(b, (const mp_limb_t*)x.value, (const mp_limb_t*)t, (const mp_limb_t*)Fp::modulus, tp);
    mpn_copyi(q, (const mp_limb_t*)b, SIZE);
    q[0]--;
    if (mpn_zero_p((const mp_limb_t*)q, SIZE) == 0) { // b != 1
        mpn_zero((mp_limb_t *)r.value, SIZE);
        // エラー投げたほうがいいかも
        return false;
    }

    bcnt = mpn_scan1(t, 0);
    mpn_rshift(q, (const mp_limb_t*)t, SIZE, bcnt);

    mpn_zero(z, SIZE);
    z[0] = 2;
    while(1) { // find quadratic non-residue
        powMod(b, (const mp_limb_t*)z, (const mp_limb_t*)t, (const mp_limb_t*)Fp::modulus, tp);
 
        b[0]++; // b + 1 == p
        if (mpn_cmp((const mp_limb_t*)b, (const mp_limb_t*)Fp::modulus, SIZE) == 0) break; // b != 1
        z[0]++;
    }

    
    powMod(c, (const mp_limb_t*)z, (const mp_limb_t*)q, (const mp_limb_t*)Fp::modulus, tp);
    powMod(t, (const mp_limb_t*)x.value, (const mp_limb_t*)q, (const mp_limb_t*)Fp::modulus, tp);

    mpn_add_1(b, (const mp_limb_t*)q, SIZE, 1);
    mpn_rshift(b, (const mp_limb_t*)b, SIZE, 1);

    powMod((mp_limb_t *)r.value, (const mp_limb_t*)x.value, (const mp_limb_t*)b, (const mp_limb_t*)Fp::modulus, tp);
    
    while(1) {
        mpn_copyi(z, (const mp_limb_t*)t, SIZE);
        z[0]--;
        if (mpn_zero_p(t, SIZE) == 1) { // t == 0
            mpn_zero((mp_limb_t *)r.value, SIZE);
            return false;
        } else if (mpn_zero_p(z, SIZE) == 1) { // t == 1
            if ((r.value[0] & 1) == 1) {
                mpn_sub_n((mp_limb_t *)r.value, (const mp_limb_t *)Fp::modulus,  
                        (const mp_limb_t*)r.value, SIZE);
            }
            return true;
        }
        mpn_copyi(z, (const mp_limb_t *)t, SIZE);
        unsigned int i = 1;
        mp_limb_t tmp[SIZE*2];
        mp_limb_t tmp_q[SIZE+1];

        mpn_copyi(b, (const mp_limb_t*)z, SIZE);
        b[0]--;
        while(mpn_zero_p(b, SIZE) == 0) {
            sqrMod(z, (const mp_limb_t*)z, (const mp_limb_t*)Fp::modulus, tmp, tmp_q);
            mpn_copyi(b, (const mp_limb_t*)z, SIZE);
            b[0]--;
            i++;
        }
        mpn_copyi(b, (const mp_limb_t *)c, SIZE);
        for(unsigned int j = 0; j < bcnt-i-1; j++) {
            sqrMod(b, (const mp_limb_t*)b, (const mp_limb_t*)Fp::modulus, tmp, tmp_q);
        }
        sqrMod(c, (const mp_limb_t*)c, (const mp_limb_t*)Fp::modulus, tmp, tmp_q);
        mulMod(t, (const mp_limb_t*)t, (const mp_limb_t*)c, (const mp_limb_t*)Fp::modulus, tmp, tmp_q);
        mulMod((mp_limb_t*)r.value, (const mp_limb_t*)r.value, (const mp_limb_t*)b, 
                (const mp_limb_t*)Fp::modulus, tmp, tmp_q);
    }
}

static inline void powMod(mp_limb_t* r, const mp_limb_t* x, const mp_limb_t* e, const mp_limb_t* modulus,
        mp_limb_t* tp) {
    mpn_sec_powm(r, x, SIZE, e, SIZE*GMP_NUMB_BITS, modulus, SIZE, tp);
}

static inline void sqrMod(mp_limb_t* r, const mp_limb_t* x, const mp_limb_t* modulus, 
        mp_limb_t* tmp,  mp_limb_t* q) {
    mpn_sqr(tmp, x, SIZE);
    mpn_tdiv_qr(q, r, 0, (const mp_limb_t*)tmp, SIZE*2, modulus, SIZE);
}
static inline void mulMod(mp_limb_t* z, const mp_limb_t* x, const mp_limb_t* y, const mp_limb_t* modulus, 
        mp_limb_t* tmp,  mp_limb_t* q) {
    mpn_mul_n(tmp, x, y, SIZE);
    mpn_tdiv_qr(q, z, 0, (const mp_limb_t*)tmp, SIZE*2, modulus, SIZE);
} 

#ifdef SECP521
static inline void mod(mp_limb_t *z, const mp_limb_t *XY, const mp_limb_t *p, mp_limb_t *t, mp_limb_t *s) {
    // (T + (T mod R)*N) / R
    mpn_zero(t, SIZE*2);
    mpn_zero(s, SIZE*2);
    mpn_and_n(t, XY, p, SIZE); // T mod R
    for (size_t i = 0; i < SIZE; i++) {
        s[i+8] = t[i];
    }
    mpn_lshift(s, (const mp_limb_t*)s, SIZE*2, 9);
    sub_n(s, s, t, SIZE*2);
    add_n(t, s, (mp_limb_t *)XY, SIZE*2); // (T + (T mod R)*N)

    mpn_rshift(t, (const mp_limb_t*)t, SIZE*2, 9);
    for (size_t i = 0; i < SIZE; i++) { // (T + (T mod R)*N) / R
        t[i] = t[i+8];
    }
    if (mpn_cmp((const mp_limb_t*)t, p, SIZE) >= 0) {
        sub_n(t, t, (mp_limb_t*)p, SIZE);
    }
    mpn_copyi(z, (const mp_limb_t*)t, SIZE);
}
#endif
