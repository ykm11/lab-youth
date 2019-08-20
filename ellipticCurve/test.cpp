#include<iostream>
//#include "curve.h"
#include<time.h>

#include "FP.h"
//#include<gmpxx.h>

#if 1
void benchmark() {
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);
    std::cout << Fp::modulus << "\n";
    Fp x, y, z;

    x = Fp(mpz_class("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10)); 
    y = Fp(mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16));

    const int n = 100000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        add(x, x, y);
    }
    time_t end = clock();
    printf("time = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}
#endif

/*
void order_test() {
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class n = mpz_class("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10); 
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC.point(gx, gy);
    Point R;
    R = G*n;

    if (R.x == zero && R.y == one && R.z == zero) {
        std::cout << "order test: OK" << std::endl;
    } else {
        std::cout << "order test: FAILED" << std::endl;
    }
}

void isEqual_test() {
    Fp::setModulo(19);

    Fp x = Fp(12);
    Fp y = Fp(12 - 19);
    Fp z = Fp(-12);

    std::cout << "Fp isEq test1: ";
    if (x == y) {
        std::cout << "OK" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
    }

    std::cout << "Fp isEq test2: ";
    if (x != z) {
        std::cout << "OK" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
    }

}

*/
int main() {
    //order_test();
    //isEqual_test();
    benchmark();
}
