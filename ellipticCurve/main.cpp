#include<iostream>
#include "curve.h"


int main() {
    mpz_class p = 65537;
    Fp::setModulo(p);

    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("0", 10);

    EllipticCurve EC = EllipticCurve(a, b);

    Fp x, y, inv_y, inv_g;
    Fp g, d;

    Point P = EC(2, 8160);
    print(P);
    P.xy(x, y);
    invmod(inv_y, y); // inv_y := y^{-1}
    mul(g, x, inv_y); // x / y

    Point R = P*18;
    print(R);
    R.xy(x, y);
    invmod(inv_y, y); 
    mul(y, x, inv_y); // x / y

    invmod(inv_g, g); 
    mul(d, y, inv_g);
    std::cout << d.value << std::endl;
}
