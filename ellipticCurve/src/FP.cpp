#include "FP.h"
#include <gmpxx.h>
#include <iostream>

#ifndef USE_MPN
mpz_class Fp::modulus;
#else
mp_limb_t Fp::modulus[SIZE];
#endif

void Fp::mulInt(Fp& z, const Fp& x, int scalar) {
#ifndef USE_MPN
    z.value = x.value * scalar;
    z.value %= modulus;
#else
    mp_limb_t tmp_z[SIZE + 1] = {0};
    mp_limb_t q[2] = {0};
    mpn_copyi((mp_limb_t *)tmp_z, (const mp_limb_t *)x.value, SIZE);

    mpn_mul_1(tmp_z, (const mp_limb_t *)tmp_z, SIZE + 1, (mp_limb_t)scalar);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, SIZE + 1, (const mp_limb_t *)Fp::modulus, SIZE);
#endif
}

void Fp::setModulo(const mpz_class& v) {
#ifndef USE_MPN
    modulus = v;
#else
    getArray(modulus, SIZE, v, v.get_mpz_t()->_mp_size);
#endif
}

#ifdef USE_MPN
void Fp::setModulo(const uint64_t p[SIZE]) {
    mpn_copyi((mp_limb_t *)modulus, (const mp_limb_t *)p, SIZE);
}
#endif


void Fp::neg(Fp& r, const Fp& x) {
#ifndef USE_MPN
    mpz_sub(r.value.get_mpz_t(), Fp::modulus.get_mpz_t(), x.value.get_mpz_t());
#else
    mpn_sub_n((mp_limb_t *)r.value, (const mp_limb_t *)Fp::modulus, (const mp_limb_t *)x.value, SIZE);
#endif
}

bool isEq(const Fp& x, const Fp& y) {
#ifndef USE_MPN
    return x.value == y.value;
#else
    return mpn_cmp((const mp_limb_t*)x.value, (const mp_limb_t*)y.value, SIZE) == 0;
#endif
}

void add(Fp& z, const Fp& x, const Fp& y) {
#ifndef USE_MPN
    z.value = x.value + y.value;
    if(z.value >= Fp::modulus) {
        z.value -= Fp::modulus;
    }
#else
    if (mpn_add_n((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE)) {
        mp_limb_t r[SIZE] = {0};
        mpn_sub_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
        mpn_copyi((mp_limb_t *)z.value, (const mp_limb_t *)r, SIZE);
        return;
    }

    if (mpn_cmp((const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE) >= 0) {
        mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
    }
#endif
}

#ifdef USE_MPN
void add(Fp& z, const Fp& x, uint64_t scalar) {
    if (mpn_add_1((mp_limb_t *)z.value, (const mp_limb_t *)x.value, SIZE, (mp_limb_t)scalar)) {
        mp_limb_t r[SIZE] = {0};
        mpn_sub_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
        mpn_copyi((mp_limb_t *)z.value, (const mp_limb_t *)r, SIZE);
        return;
    }
    if (mpn_cmp((const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE) >= 0) {
        mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
    }
}
#endif

void sub(Fp& z, const Fp& x, const Fp& y) {
#ifndef USE_MPN
    z.value = x.value - y.value;
    if (z.value < 0) {
        z.value += Fp::modulus;
    }
#else
    if (mpn_sub_n((mp_limb_t *)z.value, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE)) {
        mp_limb_t r[SIZE] = {0};
        mpn_add_n(r, (const mp_limb_t *)z.value, (const mp_limb_t *)Fp::modulus, SIZE);
        mpn_copyi((mp_limb_t *)z.value, r, SIZE);
    }
#endif
}

void mul(Fp& z, const Fp& x, const Fp& y) {
#ifdef SECP521
    mp_limb_t tmp_z[SIZE * 2] = {0};
    mp_limb_t t[SIZE*2] = {0};
    mp_limb_t s[SIZE*2] = {0};

    mpn_mul_n(tmp_z, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE);
    mod((mp_limb_t*)z.value, (const mp_limb_t *)tmp_z, (const mp_limb_t *)Fp::modulus, t, s);

#elif defined(USE_MPN)

    mp_limb_t tmp_z[SIZE * 2] = {0};
    mp_limb_t q[SIZE + 1] = {0};

    mpn_mul_n(tmp_z, (const mp_limb_t *)x.value, (const mp_limb_t *)y.value, SIZE);
    mpn_tdiv_qr(q, (mp_limb_t *)z.value, 0,
            tmp_z, SIZE*2, (const mp_limb_t *)Fp::modulus, SIZE);
#else
    z.value = (x.value * y.value) % Fp::modulus;
#endif
}


void invmod(Fp& r, const Fp& x) {
#ifndef USE_MPN
    mpz_invert(r.value.get_mpz_t(), x.value.get_mpz_t(), Fp::modulus.get_mpz_t());
#else
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
#endif
}

void sqr(Fp &r, const Fp &x) { // r <- x^2
#ifdef SECP521
    mp_limb_t tmp_r[SIZE * 2] = {0};
    mp_limb_t t[SIZE*2];
    mp_limb_t s[SIZE*2];

    mpn_sqr(tmp_r, (const mp_limb_t *)x.value, SIZE);
    mod((mp_limb_t*)r.value, (const mp_limb_t *)tmp_r, (const mp_limb_t *)Fp::modulus, t, s);
#elif defined(USE_MPN)
    mp_limb_t tmp_r[SIZE * 2] = {0};
    mp_limb_t q[SIZE + 1] = {0};

    mpn_sqr(tmp_r, (const mp_limb_t *)x.value, SIZE);
    mpn_tdiv_qr(q, (mp_limb_t *)r.value, 0,
            tmp_r, SIZE*2, (const mp_limb_t *)Fp::modulus, SIZE);
#else
    mpz_powm_ui(r.value.get_mpz_t(), x.value.get_mpz_t(), 2, Fp::modulus.get_mpz_t());
#endif
}

#ifndef USE_MPN
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
#endif
