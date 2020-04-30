#pragma once
#include "FP.h"
#include <gmpxx.h>
#include <iostream>

mp_limb_t Fp::modulus[YKM_ECC_MAX_SIZE];

void Fp::mulInt(Fp& z, const Fp& x, int scalar) {
    mp_limb_t tmp_z[size + 1];
    mp_limb_t q[2];
    
    mpn_zero(tmp_z, size+1);
    mpn_copyi((mp_limb_t *)tmp_z, (const mp_limb_t *)x.value, size);

    mpn_mul_1(tmp_z, (const mp_limb_t *)tmp_z, size + 1, (mp_limb_t)scalar);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, size + 1, (const mp_limb_t *)Fp::modulus, size);
}

void Fp::setModulo(const mpz_class& v) {
    getArray(modulus, size, v, v.get_mpz_t()->_mp_size);
}

void Fp::setModulo(const mp_limb_t p[YKM_ECC_MAX_SIZE]) {
    mpn_copyi((mp_limb_t *)modulus, (const mp_limb_t *)p, size);
}


void Fp::neg(Fp& r, const Fp& x) {
    sub_n(r.value, Fp::modulus, (mp_limb_t*)x.value, size);
    //mpn_sub_n((mp_limb_t *)r.value, (const mp_limb_t *)Fp::modulus, (const mp_limb_t *)x.value, size);
}

bool isEq(const Fp& x, const Fp& y) {
    return mpn_cmp((const mp_limb_t*)x.value, (const mp_limb_t*)y.value, Fp::size) == 0;
}

void add(Fp& z, const Fp& x, const Fp& y) {
    if (mpn_add_n((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, Fp::size)) {
        mp_limb_t r[Fp::size];
        //mpn_sub_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, Fp::size);
        sub_n(r, z.value, Fp::modulus, Fp::size);
        mpn_copyi((mp_limb_t *)z.value, (const mp_limb_t *)r, Fp::size);
        return;
    }

    if (mpn_cmp((const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, Fp::size) >= 0) {
        //mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, Fp::size);
        sub_n(z.value, z.value, Fp::modulus, Fp::size);
    }
}

void add(Fp& z, const Fp& x, uint64_t scalar) {
    if (mpn_add_1((mp_limb_t *)z.value, (const mp_limb_t *)x.value, Fp::size, (mp_limb_t)scalar)) {
        mp_limb_t r[Fp::size];
        //mpn_sub_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, Fp::size);
        sub_n(r, z.value, Fp::modulus, Fp::size);
        mpn_copyi((mp_limb_t *)z.value, (const mp_limb_t *)r, Fp::size);
        return;
    }
    if (mpn_cmp((const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, Fp::size) >= 0) {
        //mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, Fp::size);
        mpn_sub_n(z.value, z.value, Fp::modulus, Fp::size);
    }
}

void sub(Fp& z, const Fp& x, const Fp& y) {
    if (mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, Fp::size)) {
        mp_limb_t r[Fp::size];
        //mpn_add_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, Fp::size);
        add_n(r, z.value, Fp::modulus, Fp::size);
        mpn_copyi((mp_limb_t *)z.value, r, Fp::size);
    }
}

void mul(Fp& z, const Fp& x, const Fp& y) {
#ifdef SECP521
    mp_limb_t tmp_z[Fp::size * 2];
    mp_limb_t t[Fp::size*2];
    mp_limb_t s[Fp::size*2];

    //mpn_mul_n(tmp_z, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, Fp::size);
    mul_n(tmp_z, (mp_limb_t *)x.value, (mp_limb_t *)y.value, Fp::size);
    mod((mp_limb_t*)z.value, (const mp_limb_t *)tmp_z, (const mp_limb_t *)Fp::modulus, t, s);

#elif defined(USE_MPN)

    mp_limb_t tmp_z[Fp::size * 2];
    mp_limb_t q[Fp::size + 1];

    //mpn_mul_n(tmp_z, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, Fp::size);
    mul_n(tmp_z, (mp_limb_t *)x.value, (mp_limb_t *)y.value, Fp::size);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, Fp::size*2, (const mp_limb_t *)Fp::modulus, Fp::size);
#else
    z.value = (x.value * y.value) % Fp::modulus;
#endif
}


void invmod(Fp& r, const Fp& x) {
    mp_limb_t g[Fp::size];
    mp_limb_t u[Fp::size];
    mp_limb_t v[Fp::size];
    mp_limb_t s[Fp::size+1];
    mp_size_t sn;

    mpn_copyi(u, (const mp_limb_t*)x.value, Fp::size);
    mpn_copyi(v, (const mp_limb_t*)Fp::modulus, Fp::size);
    mpn_gcdext(g, s, &sn, u, Fp::size, v, Fp::size);
    mpn_copyi((mp_limb_t*)r.value, (const mp_limb_t*)s, Fp::size);

    if (sn < 0) {
        //mpn_sub_n((mp_limb_t*)r.value, (const mp_limb_t*)Fp::modulus, (const mp_limb_t*)r.value, Fp::size);
        sub_n(r.value, Fp::modulus, r.value, Fp::size);
    }
}

void sqr(Fp &r, const Fp &x) { // r <- x^2
#ifdef SECP521
    mp_limb_t tmp_r[Fp::size * 2];
    mp_limb_t t[Fp::size*2];
    mp_limb_t s[Fp::size*2];

    mpn_sqr(tmp_r, (const mp_limb_t *)x.value, Fp::size);
    mod((mp_limb_t*)r.value, (const mp_limb_t *)tmp_r, (const mp_limb_t *)Fp::modulus, t, s);
#elif defined(USE_MPN)
    mp_limb_t tmp_r[Fp::size * 2];
    mp_limb_t q[Fp::size + 1];

    mpn_sqr(tmp_r, (const mp_limb_t *)x.value, Fp::size);
    mpn_tdiv_qr(q, (mp_limb_t *)r.value, 0,
            tmp_r, Fp::size*2, (const mp_limb_t *)Fp::modulus, Fp::size);
#else
    mpz_powm_ui(r.value.get_mpz_t(), x.value.get_mpz_t(), 2, Fp::modulus.get_mpz_t());
#endif
}

bool Fp::squareRoot(Fp& r, const Fp& x) {
    mp_limb_t q[size];
    mp_limb_t c[size], t[size], b[size], z[size];
    mp_limb_t tp[mpn_sec_powm_itch(size, Fp::size*GMP_NUMB_BITS, size)];
    mp_bitcnt_t bcnt;

    mpn_copyi(t, (const mp_limb_t*)Fp::modulus, size); // t = p
    mpn_rshift(t, (const mp_limb_t*)t, size, 1); // t = (p-1)/2
    powMod(b, (const mp_limb_t*)x.value, (const mp_limb_t*)t, (const mp_limb_t*)Fp::modulus, tp, size);
    mpn_copyi(q, (const mp_limb_t*)b, size);
    q[0]--;
    if (mpn_zero_p((const mp_limb_t*)q, size) == 0) { // b != 1
        mpn_zero((mp_limb_t *)r.value, size);
        // エラー投げたほうがいいかも
        return false;
    }

    bcnt = mpn_scan1(t, 0);
    mpn_rshift(q, (const mp_limb_t*)t, size, bcnt);

    mpn_zero(z, size);
    z[0] = 2;
    while(1) { // find quadratic non-residue
        powMod(b, (const mp_limb_t*)z, (const mp_limb_t*)t, (const mp_limb_t*)Fp::modulus, tp, size);
 
        b[0]++; // b + 1 == p
        if (mpn_cmp((const mp_limb_t*)b, (const mp_limb_t*)Fp::modulus, size) == 0) break; // b != 1
        z[0]++;
    }

    
    powMod(c, (const mp_limb_t*)z, (const mp_limb_t*)q, (const mp_limb_t*)Fp::modulus, tp, size);
    powMod(t, (const mp_limb_t*)x.value, (const mp_limb_t*)q, (const mp_limb_t*)Fp::modulus, tp, size);

    mpn_add_1(b, (const mp_limb_t*)q, size, 1);
    mpn_rshift(b, (const mp_limb_t*)b, size, 1);

    powMod((mp_limb_t *)r.value, (const mp_limb_t*)x.value, (const mp_limb_t*)b, 
            (const mp_limb_t*)Fp::modulus, tp, size);
    
    while(1) {
        mpn_copyi(z, (const mp_limb_t*)t, size);
        z[0]--;
        if (mpn_zero_p(t, Fp::size) == 1) { // t == 0
            mpn_zero((mp_limb_t *)r.value, size);
            return false;
        } else if (mpn_zero_p(z, size) == 1) { // t == 1
            if ((r.value[0] & 1) == 1) {
                //mpn_sub_n((mp_limb_t *)r.value, (const mp_limb_t *)Fp::modulus, (const mp_limb_t*)r.value, size);
                sub_n(r.value, Fp::modulus, r.value, size);
            }
            return true;
        }
        mpn_copyi(z, (const mp_limb_t *)t, size);
        unsigned int i = 1;
        mp_limb_t tmp[Fp::size*2];
        mp_limb_t tmp_q[Fp::size+1];

        mpn_copyi(b, (const mp_limb_t*)z, size);
        b[0]--;
        while(mpn_zero_p(b, Fp::size) == 0) {
            sqrMod(z, (const mp_limb_t*)z, (const mp_limb_t*)Fp::modulus, tmp, tmp_q, size);
            mpn_copyi(b, (const mp_limb_t*)z, size);
            b[0]--;
            i++;
        }
        mpn_copyi(b, (const mp_limb_t *)c, size);
        for(unsigned int j = 0; j < bcnt-i-1; j++) {
            sqrMod(b, (const mp_limb_t*)b, (const mp_limb_t*)Fp::modulus, tmp, tmp_q, size);
        }
        sqrMod(c, (const mp_limb_t*)c, (const mp_limb_t*)Fp::modulus, tmp, tmp_q, size);
        mulMod(t, (const mp_limb_t*)t, (const mp_limb_t*)c, (const mp_limb_t*)Fp::modulus, tmp, tmp_q, size);
        mulMod((mp_limb_t*)r.value, (const mp_limb_t*)r.value, (const mp_limb_t*)b, 
                (const mp_limb_t*)Fp::modulus, tmp, tmp_q, size);
    }
}

