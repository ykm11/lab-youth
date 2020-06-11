#pragma once
#include "FP.h"
#include <gmpxx.h>
#include <iostream>

mp_limb_t Fp::modulus[YKM_ECC_MAX_SIZE];

void Fp::mulInt(Fp& z, const Fp& x, mp_limb_t scalar) {
    mp_limb_t tmp_z[size_ + 1];
    mp_limb_t q[2];
    mp_limb_t tmp; 

    tmp = mpn_mul_1(tmp_z, x.value, size_, scalar);
    tmp_z[size_] = tmp;
    mpn_tdiv_qr(q, z.value, 0, tmp_z, size_ + 1, 
            (const mp_limb_t *)modulus, size_);
}

void Fp::setModulo(const mpz_class& v) {
    getArray(modulus, size_, v, v.get_mpz_t()->_mp_size);
}

void Fp::setModulo(const mp_limb_t p[YKM_ECC_MAX_SIZE]) {
    copy_n(modulus, (mp_limb_t *)p, size_);
}


void Fp::neg(Fp& r, const Fp& x) {
    if (zeroCmp(x)) {
        mpn_zero(r.value, Fp::size_);
        return;
    }
    sub_n(r.value, modulus, (mp_limb_t*)x.value, size_);
}

bool isEq(const Fp& x, const Fp& y) {
    return cmp_n((mp_limb_t*)x.value, (mp_limb_t*)y.value, Fp::size_) == 0;
}

void add(Fp& z, const Fp& x, const Fp& y) {
    if (add_n(z.value, (mp_limb_t *)x.value, (mp_limb_t *)y.value, Fp::size_)) {
        sub_n(z.value, z.value, Fp::modulus, Fp::size_);
        return;
    }

    if (cmp_n(z.value, Fp::modulus, Fp::size_) >= 0) {
        sub_n(z.value, z.value, Fp::modulus, Fp::size_);
    }
}

void add(Fp& z, const Fp& x, mp_limb_t scalar) {
    if (add_1(z.value, (mp_limb_t *)x.value, Fp::size_, scalar)) {
        sub_n(z.value, z.value, Fp::modulus, Fp::size_);
        return;
    }
    if (cmp_n(z.value, Fp::modulus, Fp::size_) >= 0) {
        sub_n(z.value, z.value, Fp::modulus, Fp::size_);
    }
}

void sub(Fp& z, const Fp& x, const Fp& y) {
    if (sub_n(z.value, (mp_limb_t *)x.value, (mp_limb_t *)y.value, Fp::size_)) {
        add_n(z.value, z.value, Fp::modulus, Fp::size_);
    }
}

void mul(Fp& z, const Fp& x, const Fp& y) {
#ifdef YKM_ECC_SECP521
    mp_limb_t tmp_z[Fp::size_ * 2];
    mp_limb_t t[Fp::size_*2];
    mp_limb_t s[Fp::size_*2];

    mul_n(tmp_z, (mp_limb_t *)x.value, (mp_limb_t *)y.value, Fp::size_);
    secp521Mod(z.value, tmp_z, Fp::modulus, t, s);

#elif defined(YKM_ECC_USE_MPN)
    mp_limb_t tmp_z[Fp::size_ * 2];
    mp_limb_t q[Fp::size_ + 1];

    mulMod(z.value, (mp_limb_t *)x.value, (mp_limb_t *)y.value, 
            Fp::modulus, tmp_z, q, Fp::size_);
#else
    z.value = (x.value * y.value) % Fp::modulus;
#endif
}


void invmod(Fp& r, const Fp& x) {
    mp_limb_t g[Fp::size_];
    mp_limb_t u[Fp::size_];
    mp_limb_t v[Fp::size_];
    mp_limb_t s[Fp::size_+1];
    mp_size_t sn;

    mp_size_t x_size = Fp::size_;
    while (x_size > 0) {
        if(!x.value[x_size - 1]) break;
        x_size--;
    }

    copy_n(v, (mp_limb_t *)x.value, Fp::size_);
    copy_n(u, Fp::modulus, Fp::size_);
    mpn_gcdext(g, s, &sn, u, Fp::size_, v, x_size); 
    
    mpn_zero(r.value, Fp::size_);
    if (sn == 0) { // x == 1;
        r.value[0] = 1;
        return;
    }

    mp_limb_t t[Fp::size_ * 2];
    // 1 - s*p = inv_x*x
    if (sn < 0) { // -sp > 0
        mpn_mul(t, (const mp_limb_t*)Fp::modulus, Fp::size_, (const mp_limb_t*)s, -sn); // s * p
        mpn_add_1(t, t, -sn + Fp::size_, 1); // 1 + sp
        mpn_tdiv_qr(r.value, u, 0, 
                (const mp_limb_t*)t, -sn + Fp::size_, 
                (const mp_limb_t*)x.value, x_size);

    } else { // -sp < 0
        mpn_mul(t, (const mp_limb_t*)Fp::modulus, Fp::size_, (const mp_limb_t*)s, sn); // s * p
        mpn_sub_1(t, t, sn + Fp::size_, 1); // sp - 1
        mpn_tdiv_qr(r.value, u, 0, 
                (const mp_limb_t*)t, sn + Fp::size_, 
                (const mp_limb_t*)x.value, x_size);
        sub_n(r.value, Fp::modulus, r.value, Fp::size_);
    }
}

void sqr(Fp &r, const Fp &x) { // r <- x^2
#ifdef YKM_ECC_SECP521
    mp_limb_t tmp_r[Fp::size_ * 2];
    mp_limb_t t[Fp::size_*2];
    mp_limb_t s[Fp::size_*2];

    mpn_sqr(tmp_r, (mp_limb_t *)x.value, Fp::size_);
    secp521Mod(r.value, tmp_r, Fp::modulus, t, s);
#elif defined(YKM_ECC_USE_MPN)
    mp_limb_t tmp_r[Fp::size_ * 2];
    mp_limb_t q[Fp::size_ + 1];

    sqrMod(r.value, (mp_limb_t *)x.value, Fp::modulus, tmp_r, q, Fp::size_);
#endif
}

bool Fp::squareRoot(Fp& r, const Fp& x) {
    mp_limb_t q[size_];
    mp_limb_t c[size_], t[size_], b[size_], z[size_];
    mp_limb_t tp[mpn_sec_powm_itch(size_, size_*GMP_NUMB_BITS, size_)];
    mp_bitcnt_t bcnt;

    copy_n(t, modulus, size_); // t = p
    mpn_rshift(t, (const mp_limb_t*)t, size_, 1); // t = (p-1)/2
    powMod(b, (mp_limb_t*)x.value, t, modulus, tp, size_);
    copy_n(q, b, size_);
    q[0]--;
    if (mpn_zero_p((const mp_limb_t*)q, size_) == 0) { // b != 1
        mpn_zero((mp_limb_t *)r.value, size_);
        // エラー投げたほうがいいかも
        return false;
    }

    bcnt = mpn_scan1(t, 0);
    mpn_rshift(q, (const mp_limb_t*)t, size_, bcnt);

    mpn_zero(z, size_);
    z[0] = 2;
    while(1) { // find quadratic non-residue
        powMod(b, z, t, modulus, tp, size_);
 
        b[0]++; // b + 1 == p
        if (cmp_n(b, Fp::modulus, size_) == 0) break; // b != 1
        z[0]++;
    }

    
    powMod(c, z, q, modulus, tp, size_);
    powMod(t, (mp_limb_t*)x.value, q, modulus, tp, size_);

    mpn_add_1(b, (const mp_limb_t*)q, size_, 1);
    mpn_rshift(b, (const mp_limb_t*)b, size_, 1);

    powMod(r.value, (mp_limb_t*)x.value, b, modulus, tp, size_);
    
    while(1) {
        copy_n(z, t, size_);
        z[0]--;
        if (mpn_zero_p(t, size_) == 1) { // t == 0
            mpn_zero((mp_limb_t *)r.value, size_);
            return false;
        } else if (mpn_zero_p(z, size_) == 1) { // t == 1
            if ((r.value[0] & 1) == 1) {
                sub_n(r.value, modulus, r.value, size_);
            }
            return true;
        }
        copy_n(z, t, size_);
        unsigned int i = 1;
        mp_limb_t tmp[size_*2];
        mp_limb_t tmp_q[size_+1];

        copy_n(b, z, size_);

        b[0]--;
        while(mpn_zero_p(b, size_) == 0) {
            sqrMod(z, z, modulus, tmp, tmp_q, size_);
            copy_n(b, z, size_);
            b[0]--;
            i++;
        }
        copy_n(b, c, size_);
        for(unsigned int j = 0; j < bcnt-i-1; j++) {
            sqrMod(b, b, modulus, tmp, tmp_q, size_);
        }
        sqrMod(c, c, modulus, tmp, tmp_q, size_);
        mulMod(t, t, c, modulus, tmp, tmp_q, size_);
        mulMod(r.value, r.value, b, modulus, tmp, tmp_q, size_);
    }
}

