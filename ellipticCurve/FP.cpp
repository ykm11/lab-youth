#include "FP.h"
#include <gmpxx.h>

mpz_class Fp::modulus;

void Fp::mulInt(Fp& z, const Fp& x, int scalar) {
    z.value = x.value * scalar;
    z.value %= modulus;
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
