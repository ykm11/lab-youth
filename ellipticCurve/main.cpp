#include<iostream>
#include "curve.h"

int main() {
    mpz_class a = 2;
    mpz_class b = 3;
    mpz_class p = 19;

    EllipticCurve EC = EllipticCurve(a, b, p);
    EC.printCurve();

    Point P = EC(1, 5);
    Point Q = EC(3, 6);

    Point R;
    R = P + Q;
    std::cout << R.x << " : " << R.y << " : " << R.z << std::endl;

    Fp::setModulo(23);

    Fp x = Fp(21);
    Fp y = Fp(9);

    Fp z;
    z = x + y;
    std::cout << z.value << std::endl;

    z = x - y;
    std::cout << z.value << std::endl;

    z = x * y;
    std::cout << z.value << std::endl;
}
