#ifdef USE_MIMALLOC
/*
    g++ -O3 -DNDEBUG -I /usr/local/lib/mimalloc/include /usr/local/lib/mimalloc/build/libmimalloc.a -lpthread -lgmpxx -lgmp ./main.cpp
*/
#include <gmp.h>
#include <mimalloc.h>
static struct UseMiMalloc {
    static void* mi_realloc_wrapper(void *p, size_t, size_t n)
    {
        return mi_realloc(p, n);
    }
    static void mi_free_wrapper(void *p, size_t)
    {
        mi_free(p);
    }
    UseMiMalloc()
    {
        puts("set GMP memory functions before using mpz_class");
        mp_set_memory_functions(mi_malloc, mi_realloc_wrapper, mi_free_wrapper);
    }
} g_UseMiMalloc;
#endif

#include<iostream>
#include "attack_curve.h"

mpz_class Fp::modulus;
Fp EllipticCurve::a;
Fp EllipticCurve::b;

int main() {
    mpz_class p = mpz_class("115792089237316195423570985008687907853269984665640564039457584007908834359199", 10);
    Fp::setModulo(p);

    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);

    EllipticCurve E = EllipticCurve(a, b);

    mpz_class px, py;
    mpz_class qx, qy;

    mpz_class order = mpz_class("115792089237316195423570985008687907853269984665640564039457584007908834359200", 10);

    mpz_class p_card = mpz_class("237140208281", 10);

    px = mpz_class("55066263022277343669578718895168534326250603453777594175500187360389116729240", 10);
    py = mpz_class("82420052799988522717532479648954225145345455265426509542622303173620650331589", 10);
    qx = mpz_class("56df2adff3c3749cc4c62c9e7da339dc02d157868a1d76f9d058d634d6a9525f", 16);
    qy = mpz_class("c167d7eb600437e2d6ead69ebcf2b1b2f88939c0fafda0b19aa3db33f5024b43", 16);

    Point P = E(px, py);
    Point Q = E(qx, qy);

    if ( order % p_card != 0) {
        return 1;
    }
    mpz_class beta = order / p_card;

    Point P_ = P * beta;
    Point Q_ = Q * beta;

    mpz_class d = pollard_rho_ECDLP(P_, Q_, E, p_card);
    std::cout << d << std::endl;
    std::cout << (Q_ == P_*d) << std::endl;
}
