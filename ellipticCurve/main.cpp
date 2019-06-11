#include<iostream>
#include "curve.h"

int main() {
    mpz_class a = 0;
    mpz_class b = 1;
    mpz_class p = 17;

    EllipticCurve EC = NewCurve(a, b, p);
    EC.printCurve();

    mpz_class x, y;
    x = 2;
    y = 3;
    point P = EC.Point(x, y);
    std::cout << P.x << std::endl;
    std::cout << P.y << std::endl;

    return 0;
}
