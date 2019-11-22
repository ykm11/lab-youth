#pragma once
#include<gmpxx.h>


class RSA {
public:
    mpz_class n, e;

    RSA() { }
    RSA(const mpz_class &p1, const mpz_class &p2, const mpz_class &exp) :
        n(p1*p2), e(exp), p(p1), q(p2) {
        mpz_class t;

        // GCD(e, p-1) == GCD(e, q-1) == 1のチェックは今のところ無し
        
        //t = (p - 1) * (q - 1);
        t = n - p - q + 1;
        mpz_invert(d.get_mpz_t(), e.get_mpz_t(), t.get_mpz_t()); // for (mod n)
        mpz_invert(inv_q.get_mpz_t(), q.get_mpz_t(), p.get_mpz_t());

        // for CRT-RSA
        t = p - 1;
        mpz_invert(dp.get_mpz_t(), e.get_mpz_t(), t.get_mpz_t());
        t = q - 1;
        mpz_invert(dq.get_mpz_t(), e.get_mpz_t(), t.get_mpz_t());
    }

    void encrypt(mpz_class &c, const mpz_class &m) const; // c <- m^e mod n
    void decrypt(mpz_class &m, const mpz_class &c) const; // m <- c^d mod n
    void crt_decrypt(mpz_class &m, const mpz_class &c) const; // m <- c^d mod n

private:
    mpz_class p, q, inv_q, d, dp, dq;
};
