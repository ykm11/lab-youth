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


    // Fp = 19
    // x = Fp(12)
    // y = Fp(14)

    Fp::setModulo(19);

    Fp x = Fp(-2);
    Fp y = Fp(11);

    Fp z;
    z = x + y;
    std::cout << z.value << std::endl;

    z = x - y;
    std::cout << z.value << std::endl;

    z = x * y;
    std::cout << z.value << std::endl;
}
