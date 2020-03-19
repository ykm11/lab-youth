#include "FP.h"
#include <gmpxx.h>
#include <iostream>

#ifndef USE_MPN
mpz_class Fp::modulus;

void Fp::mulInt(Fp& z, const Fp& x, int scalar) {
    z.value = x.value * scalar;
    z.value %= modulus;
}

void Fp::neg(Fp& r, const Fp& x) {
    mpz_sub(r.value.get_mpz_t(), Fp::modulus.get_mpz_t(), x.value.get_mpz_t());
}

void Fp::setModulo(const mpz_class& v) {
    modulus = v;
}

void add(Fp& z, const Fp& x, const Fp& y) {
    z.value = x.value + y.value;
    if(z.value >= Fp::modulus) {
        z.value -= Fp::modulus;
    }
}

void sub(Fp& z, const Fp& x, const Fp& y) {
    z.value = x.value - y.value;
    if (z.value < 0) {
        z.value += Fp::modulus;
    }
}

void mul(Fp& z, const Fp& x, const Fp& y) {
    z.value = (x.value * y.value) % Fp::modulus;
}

void invmod(Fp& r, const Fp& x) {
    mpz_invert(r.value.get_mpz_t(), x.value.get_mpz_t(), Fp::modulus.get_mpz_t());
}

bool isEq(const Fp& x, const Fp& y) {
    return x.value == y.value;
}

void sqr(Fp &r, const Fp &x) { // r <- x^2
    mpz_powm_ui(r.value.get_mpz_t(), x.value.get_mpz_t(), 2, Fp::modulus.get_mpz_t());
}

bool Fp::squareRoot(Fp& r, const Fp& x) {
    mpz_class q;
    mpz_class c, t, b, z;
    mp_bitcnt_t bcnt = 1;

    mpz_sub_ui(t.get_mpz_t(), Fp::modulus.get_mpz_t(), 1); // t = p - 1
    mpz_tdiv_q_2exp(t.get_mpz_t(), t.get_mpz_t(), 1); // t = (p-1)/2
    powMod(b, x.value, t, Fp::modulus);
    if (b != 1) {
        r.value = 0;
        // エラー投げたほうがいいかも
        return false;
    }

    while(mpz_divisible_2exp_p(t.get_mpz_t(), bcnt+1)) {
        bcnt++;
    }
    mpz_tdiv_q_2exp(q.get_mpz_t(), t.get_mpz_t(), bcnt-1); // p-1 == q * 2^{bcnt}

    z = 2;
    while(1) { // find quadratic non-residue
        powMod(b, z, t, Fp::modulus);
        if (b != 1) break;
        mpz_add_ui(z.get_mpz_t(), z.get_mpz_t(), 1);
    }

    powMod(c, z, q, Fp::modulus);
    powMod(t, x.value, q, Fp::modulus);
    
    mpz_add_ui(b.get_mpz_t(), q.get_mpz_t(), 1);
    mpz_tdiv_q_2exp(b.get_mpz_t(), b.get_mpz_t(), 1); // (q+1)/2
    powMod(r.value, x.value, b, Fp::modulus);
    
    while(1) {
        if (t == 0) {
            r.value = 0;
            return false;
        } else if (t == 1) {
            if (mpz_tstbit(r.value.get_mpz_t(), 0) == 1) {
                Fp::neg(r, r);
            }
            return true;
        }
        z = std::move(t);
        unsigned int i;
        for (i = 1; z != 1; i++) {
            sqrMod(z, z, Fp::modulus);
        }
        b = std::move(c);
        for(unsigned int j = 0; j < bcnt-i-1; j++) {
            sqrMod(b, b, Fp::modulus);
        }
        sqrMod(c, b, Fp::modulus);
        mulMod(t, t, c, Fp::modulus);
        mulMod(r.value, r.value, b, Fp::modulus);
    }
}

#else

#define SIZE 4
uint64_t Fp::modulus[SIZE];
void Fp::setModulo(const uint64_t p[SIZE]) {
    mpn_copyi((mp_limb_t *)modulus, (const mp_limb_t *)p, SIZE);
}

void Fp::setModulo(const mpz_class& v) {
    mpz_class x = v;
    for (size_t i = 0; i < SIZE; i++) {
        modulus[i] = mpz_get_ui(x.get_mpz_t());
        mpz_tdiv_q_2exp(x.get_mpz_t(), x.get_mpz_t(), 64);
    }
}

void add(Fp& z, const Fp &x, const Fp &y) {
    mp_limb_t carry = mpn_add_n((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE);
    if (carry == 1) {
        mp_limb_t r[SIZE+1] = {0};
        mp_limb_t p[SIZE+1] = {0};
        r[SIZE] = 1;

        mpn_copyi(r, (const mp_limb_t *)z.value, SIZE);
        mpn_copyi(p, (const mp_limb_t *)Fp::modulus, SIZE);
        mpn_sub_n(r, (const mp_limb_t *)r, (const mp_limb_t *)p, SIZE+1);
        mpn_copyi((mp_limb_t *)z.value, (const mp_limb_t *)r, SIZE);
        return;
    }

    if (mpn_cmp((const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE) >= 0) {
        mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
    }
}

void add(Fp& z, const Fp& x, uint64_t scalar) {
    mpn_add_1((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (mp_limb_t)scalar, SIZE);
    if (mpn_cmp((const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE) >= 0) {
        mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
    }
}

void sub(Fp& z, const Fp &x, const Fp &y) {
    if (mpn_cmp((const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE) < 0) { // x < y
        mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)y.value, (const mp_limb_t *)x.value, SIZE);
        mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, (const mp_limb_t *)z.value, SIZE);
        return;
    }
    mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE);
}

void mul(Fp& z, const Fp &x, const Fp &y) {
    mp_limb_t tmp_z[SIZE * 2] = {0};
    mp_limb_t q[SIZE + 1] = {0};

    mpn_mul_n(tmp_z, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, SIZE*2, (const mp_limb_t *)Fp::modulus, SIZE);
}

void sqr(Fp& z, const Fp &x) {
    mp_limb_t tmp_z[SIZE * 2] = {0};
    mp_limb_t q[SIZE + 1] = {0};

    mpn_sqr(tmp_z, (const mp_limb_t *)x.value, SIZE);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, SIZE*2, (const mp_limb_t *)Fp::modulus, SIZE);
}

void Fp::mulInt(Fp &z, const Fp &x, int scalar) {
    mp_limb_t tmp_z[SIZE + 1] = {0};
    mp_limb_t q[2] = {0};
    mpn_copyi((mp_limb_t *)tmp_z, (const mp_limb_t *)x.value, SIZE);

    mpn_mul_1(tmp_z, (const mp_limb_t *)tmp_z, SIZE + 1, (mp_limb_t)scalar);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, SIZE + 1, (const mp_limb_t *)Fp::modulus, SIZE);
}

void Fp::neg(Fp &z, const Fp &x) {
    mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, (const mp_limb_t *)x.value, SIZE);
}

bool isEq(const Fp& x, const Fp& y) {
    return mpn_cmp((const mp_limb_t*)x.value, (const mp_limb_t*)y.value, SIZE) == 0;
}

void invmod(Fp& r, const Fp& x) {
    mpz_t mr, mx, modulus;
    mpz_init(mr);
    set_mpz_t(mx, x.value, SIZE);
    set_mpz_t(modulus, Fp::modulus, SIZE);
    mpz_invert(mr, mx, modulus);

    for (size_t i = 0; i < SIZE; i++) {
        r.value[i] = mpz_get_ui(mr);
        mpz_tdiv_q_2exp(mr, mr, 64);
    }
}

#endif
