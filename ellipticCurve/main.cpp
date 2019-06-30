#include<iostream>
#include "curve.h"

int main() {
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

    z = z * two;
    std::cout << z.value << std::endl;
}
