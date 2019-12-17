#pragma once

#include<gmpxx.h>
#include<time.h>
#include<iostream>

extern mpz_class n, g;
extern mpz_class n_, r, r_; // for MR
extern mpz_class r2;

extern mpz_class X, Y, XY, S;
extern mpz_class T;

void setup();
void MontRe(mpz_class &S, const mpz_class &XY);
void MMM(mpz_class &ret, const mpz_class &x, const mpz_class &y);
void powm(mpz_class &R, const mpz_class &base, const mpz_class &exp, const mpz_class &N);
void powmMont(mpz_class &R, const mpz_class &base, const mpz_class &exp);
void powm_slide(mpz_class &R, const mpz_class &base, const mpz_class &exp, const mpz_class &N);

inline void sqrMod(mpz_class &z, const mpz_class x, const mpz_class m) {
    mpz_powm_ui(z.get_mpz_t(), x.get_mpz_t(), 2, m.get_mpz_t());
}

inline void mulMod(mpz_class &z, const mpz_class &x, const mpz_class &y, const mpz_class &m) {
    mpz_mul(z.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
    mpz_mod(z.get_mpz_t(), z.get_mpz_t(), m.get_mpz_t());
}
