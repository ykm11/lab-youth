#include<iostream>
#include "curve.h"


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

    mpz_class p_card = mpz_class("2", 10);

    px = mpz_class("55066263022277343669578718895168534326250603453777594175500187360389116729240", 10);
    py = mpz_class("82420052799988522717532479648954225145345455265426509542622303173620650331589", 10);
    qx = mpz_class("56df2adff3c3749cc4c62c9e7da339dc02d157868a1d76f9d058d634d6a9525f", 16);
    qy = mpz_class("c167d7eb600437e2d6ead69ebcf2b1b2f88939c0fafda0b19aa3db33f5024b43", 16);

    Point P = E(px, py);
    Point Q = E(qx, qy);

    mpz_class beta = order / p_card;

    Point P_ = P * beta;
    Point Q_ = Q * beta;

    mpz_class d = pollard_rho_ECDLP(P_, Q_, E, p_card);
    std::cout << d << std::endl;
    std::cout << (Q_ == P_*d) << std::endl;
}
