#include "FP.h"
#include <gmpxx.h>

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

void mul(Fp& z, const Fp& x, int scalar) {
    z.value = (x.value * scalar);
    z.value %= Fp::modulus;
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

