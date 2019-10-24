#include<gmpxx.h>
#include "rsa.h"

    
void RSA::encrypt(mpz_class &c, const mpz_class &m) const { // c <- m^e mod n
    mpz_powm(c.get_mpz_t(), m.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
}
void RSA::decrypt(mpz_class &m, const mpz_class &c) const { // m <- c^d mod n
    mpz_powm(m.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
}
void RSA::crt_decrypt(mpz_class &m, const mpz_class &c) const { // m <- c^d mod n
    mpz_class mq;

    mpz_powm(m.get_mpz_t(), c.get_mpz_t(), dp.get_mpz_t(), p.get_mpz_t()); // m <- C^{dp} mod p
    mpz_powm(mq.get_mpz_t(), c.get_mpz_t(), dq.get_mpz_t(), q.get_mpz_t()); // mq <- C^{dq} mod q

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
    m = mq + (((m - mq) * inv_q) % p) * q;
    if (m < 0) {
        m = m + n;
    }
#endif
}

