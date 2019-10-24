#pragma once
#include<gmpxx.h>

class RSA {
public:
    mpz_class n, e;

    RSA() { }
    RSA(mpz_class &p1, mpz_class &p2, mpz_class &exp) :
        n(p1*p2), e(exp), p(p1), q(p2) {
        mpz_class t;

        // GCD(e, p-1) == GCD(e, q-1) == 1のチェックは今のところ無し
        
        t = (p - 1) * (q - 1); 
        mpz_invert(d.get_mpz_t(), e.get_mpz_t(), t.get_mpz_t()); // for (mod n)
        mpz_invert(inv_q.get_mpz_t(), q.get_mpz_t(), p.get_mpz_t());

        // for CRT-RSA
        t = p - 1;
        mpz_invert(dp.get_mpz_t(), e.get_mpz_t(), t.get_mpz_t());
        t = q - 1;
        mpz_invert(dq.get_mpz_t(), e.get_mpz_t(), t.get_mpz_t());
    }

    void encrypt(mpz_class &c, const mpz_class &m) const { // c <- m^e mod n
        mpz_powm(c.get_mpz_t(), m.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
    }
    void decrypt(mpz_class &m, const mpz_class &c) const { // m <- c^d mod n
        mpz_powm(m.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
    }
    void crt_decrypt(mpz_class &m, const mpz_class &c) const { // m <- c^d mod n
        mpz_class mp, mq;

        mpz_powm(mp.get_mpz_t(), c.get_mpz_t(), dp.get_mpz_t(), p.get_mpz_t());
        mpz_powm(mq.get_mpz_t(), c.get_mpz_t(), dq.get_mpz_t(), q.get_mpz_t());
#if 0
        mpz_sub(m.get_mpz_t(), mp.get_mpz_t(), mq.get_mpz_t()); // mp - mq
        mpz_mul(m.get_mpz_t(), m.get_mpz_t(), inv_q.get_mpz_t()); // (mp - mq) * inv_q
        mpz_mod(m.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t()); // ((mp - mq) * inv_q) % p
        mpz_addmul(mq.get_mpz_t(), m.get_mpz_t(), q.get_mpz_t()); // mq + (((mp - mq) * inv_q) % p)*q
        if (mq < 0) {
            mpz_add(m.get_mpz_t(), mq.get_mpz_t(), n.get_mpz_t());
            //m = mq + n;
        } else {
            m = std::move(mq);
            //m = mq;
        }
#else
        m = mq + (((mp - mq) * inv_q) % p) * q;
        if (m < 0) {
            m = m + n;
        }
#endif
    }

private:
    mpz_class p, q, inv_q, d, dp, dq;
};
