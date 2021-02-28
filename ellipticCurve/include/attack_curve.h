#pragma once
#include "curve.h"
#include "FP.h"


// Attacks against weak elliptic curve
//
void pollard_rho_f(const Point& alpha, const Point& beta, Point& x,
        mpz_class& a, mpz_class& b, const mpz_class& order) {
    Fp x_x, x_y;
    if (!zeroCmp(x.z)) x.xy(x_x, x_y);

    if (x_x.value[0] % 3 == 0) {
        add(x, beta, x);

        b = (b + 1) % order;
    } else if (x_x.value[0] % 3 == 1) {
        EllipticCurve::dbl(x, x);
        a = (a * 2) % order;
        b = (b * 2) % order;
    } else if(x_x.value[0] % 3 == 2) {
        add(x, alpha, x);
        a = (a + 1) % order;

    }
}

mpz_class pollard_rho_ECDLP(const Point& alpha, const Point& beta,
        const EllipticCurve& curve, const mpz_class& order) {
    // beta = [d] * alpha
    mpz_class d = 0;
    mpz_class a, b, A, B;

    a = 0; b = 0;
    A = 0; B = 0;
    Point x(0,1,0);
    Point X(0,1,0);

    for(int i = 0; i < order; i++) {
        pollard_rho_f(alpha, beta, x, a, b, order);
        pollard_rho_f(alpha, beta, X, A, B, order);
        pollard_rho_f(alpha, beta, X, A, B, order);

        if (x == X && !(A == 0 && a == 0 && b == 0 && B == 0)) {
            d = (B - b) % order;
            mpz_invert(d.get_mpz_t(), d.get_mpz_t(), order.get_mpz_t());
            d = (d * (a - A)) % order;
            if (d < 0) {
                d = d + order;
            }
            return d;
        }
    }
    return d;
}
