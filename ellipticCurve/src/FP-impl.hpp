#pragma once
#include "FP.h"
#include <gmpxx.h>
#include <iostream>

mpz_class Fp::modulus;

void Fp::mulInt(Fp& z, const Fp& x, int scalar) {
    z.value = x.value * scalar;
    z.value %= modulus;
}

void Fp::setModulo(const mpz_class& v) {
    modulus = v;
}


void Fp::neg(Fp& r, const Fp& x) {
    if (x.value == 0) {
        return;
    }
    r.value = modulus - x.value;
}

bool isEq(const Fp& x, const Fp& y) {
    return x.value == y.value;
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

void sqr(Fp &r, const Fp &x) { // r <- x^2
    r.value = (x.value * x.value) % Fp::modulus;
}

bool Fp::squareRoot(Fp& r, const Fp& x) {
    mpz_class q;
    mpz_class c, t, b, z;
    mp_bitcnt_t bcnt = 1;

    t = modulus - 1; // t = p - 1
    t /= 2;
    powMod(b, x.value, t, modulus);
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
        powMod(b, z, t, modulus);
        if (b != 1) break;
        mpz_add_ui(z.get_mpz_t(), z.get_mpz_t(), 1);
    }

    powMod(c, z, q, modulus);
    powMod(t, x.value, q, modulus);
    
    b = q + 1;
    b /= 2;
    powMod(r.value, x.value, b, modulus);
    
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
            sqrMod(z, z, modulus);
        }
        b = std::move(c);
        for(unsigned int j = 0; j < bcnt-i-1; j++) {
            sqrMod(b, b, modulus);
        }
        sqrMod(c, b, modulus);
        mulMod(t, t, c, modulus);
        mulMod(r.value, r.value, b, modulus);
    }
}

