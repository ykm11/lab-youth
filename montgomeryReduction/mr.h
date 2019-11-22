#pragma once

#include<gmpxx.h>
#include<time.h>
#include<iostream>

mpz_class n, g;
mpz_class n_, r, r_; // for MR
mpz_class r2;

mpz_class X, Y, XY, S;
mpz_class T;


void setup() {
    size_t n_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    mpz_ui_pow_ui(r.get_mpz_t(), 2, n_bits); // r := n^{k+1} ( > n)
    mpz_gcdext(g.get_mpz_t(), n_.get_mpz_t(), r_.get_mpz_t(), n.get_mpz_t(), r.get_mpz_t()); // -n_*n + r_*r = 1
    mpz_mod(n_.get_mpz_t(), n_.get_mpz_t(), r.get_mpz_t()); // -n_ * n = 1 mod r
    mpz_mod(r_.get_mpz_t(), r_.get_mpz_t(), n.get_mpz_t()); // r_ * r  = 1 mod n
    if (n_ > 0) {
        mpz_sub(n_.get_mpz_t(), r.get_mpz_t(), n_.get_mpz_t()); // n_ <- -1 * n_ = r - n_
    } else {
        mpz_mul_ui(n_.get_mpz_t(), n_.get_mpz_t(), -1); // n_ <- -1 * -n_
    }
    mpz_sub_ui(T.get_mpz_t(), r.get_mpz_t(), 1); // R - 1 = 2^{k} - 1
    mpz_powm_ui(r2.get_mpz_t(), r.get_mpz_t(), 2, n.get_mpz_t()); // r2 <- r^{2} mod n
}

void MontRe(mpz_class &S, const mpz_class &XY) { // r <- XY*R^{-1}
    // S <- (XY + n * (n_ * XY mod R)) / R
    mpz_mul(S.get_mpz_t(), n_.get_mpz_t(), XY.get_mpz_t()); // n_ * XY
    mpz_and(S.get_mpz_t(), S.get_mpz_t(), T.get_mpz_t()); // mod 2^{k} = and (2^k - 1)
    mpz_mul(S.get_mpz_t(), S.get_mpz_t(), n.get_mpz_t()); // n * (n_ * XY mod R)
    mpz_add(S.get_mpz_t(), S.get_mpz_t(), XY.get_mpz_t()); // (XY + n * (n_ * XY mod R))
    mpz_tdiv_q_2exp(S.get_mpz_t(), S.get_mpz_t(), mpz_sizeinbase(n.get_mpz_t(), 2)); // (XY + n * (n_ * XY mod R)) / R

    if (S > n) {
        mpz_sub(S.get_mpz_t(), S.get_mpz_t(), n.get_mpz_t());
    }
}


void MMM(mpz_class &ret, const mpz_class &x, const mpz_class &y) { // r <- x*y mod n
#if 1
    mpz_class s;
    mpz_mul(s.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
    MontRe(ret, s);
    mpz_mul(s.get_mpz_t(), ret.get_mpz_t(), r2.get_mpz_t());
    MontRe(ret, s);
#else
    mpz_class s, t;
    mpz_mul(ret.get_mpz_t(), x.get_mpz_t(), r2.get_mpz_t());
    MontRe(s, ret); // MR(x * R2)
    mpz_mul(ret.get_mpz_t(), y.get_mpz_t(), r2.get_mpz_t());
    MontRe(t, ret); // MR(y * R2)
    mpz_mul(ret.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t()); // s * t
    MontRe(s, ret); // MR(s * t)
    MontRe(ret, s); // MR(s * t)
#endif
}


void powm(mpz_class &R, const mpz_class &base, const mpz_class &exp, const mpz_class &N) {
    // R <- base^{exp} mod N
    size_t e_bits = mpz_sizeinbase(exp.get_mpz_t(), 2);
#if 0
    mpz_class A, tmp;
    mpz_mul(tmp.get_mpz_t(), base.get_mpz_t(), r2.get_mpz_t()); // base * R2
    MontRe(A, tmp); // A := MR(base * R2)

    tmp = A;
    for(int i = e_bits - 2; i >= 0; i--) {
        mpz_mul(tmp.get_mpz_t(), tmp.get_mpz_t(), tmp.get_mpz_t());
        MontRe(R, tmp);

        tmp = std::move(R);
        if (mpz_tstbit(exp.get_mpz_t(), i) == 1) {
            mpz_mul(tmp.get_mpz_t(), tmp.get_mpz_t(), A.get_mpz_t());
            MontRe(R, tmp);

            tmp = std::move(R);
        }
    }
    MontRe(R, tmp);
#else
    R = base;
    for(int i = e_bits - 2; i >= 0; i--) {
        mpz_mul(R.get_mpz_t(), R.get_mpz_t(), R.get_mpz_t());
        mpz_mod(R.get_mpz_t(), R.get_mpz_t(), N.get_mpz_t());

        if (mpz_tstbit(exp.get_mpz_t(), i) == 1) {
            mpz_mul(R.get_mpz_t(), R.get_mpz_t(), base.get_mpz_t());
            mpz_mod(R.get_mpz_t(), R.get_mpz_t(), N.get_mpz_t());
        }
    }
#endif
}
