#include<iostream>
#include "curve.h"

int main() {
    mpz_class a = 2;
    mpz_class b = 3;
    mpz_class p = 19;

    EllipticCurve EC = EllipticCurve(a, b, p);
    EC.printCurve();

    EllipticCurve::Point P = EC(1, 5);
    EllipticCurve::Point Q = EC(3, 6);

    EllipticCurve::Point R;
    R = P + Q;
    std::cout << R.x << " : " << R.y << " : " << R.z << std::endl;
    return 0;
}
