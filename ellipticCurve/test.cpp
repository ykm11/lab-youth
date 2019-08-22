#include<iostream>
#include "curve.h"
#include<time.h>

#include "FP.h"


void benchmark_fp();
void isEqual_fp_test();
void order_test();


Fp EllipticCurve::a;
Fp EllipticCurve::b;

void benchmark_ec() {
    std::cout << "[*] EC benchmark\n";
    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class q = mpz_class("117289373161954235709850086879078528375642790749043841647", 10);
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC.point(gx, gy);
    Point R;
    add(R, G, G);
    const int n = 100000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        //mul(R, G, q);
        add(R, R, G);
    }
    time_t end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_fp() {
    std::cout << "[*] Fp benchmark\n";
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);
    Fp x, y, z;

    x = Fp("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10); 
    y = Fp("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);

    const int n = 100000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        invmod(x, Fp::modulus);
    }
    time_t end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

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

    Point O = EC.point(0, 1, 0);
    if (R == O) {
        std::cout << "[*] order test: OK" << std::endl;
    } else {
        std::cout << "[*] order test: FAILED" << std::endl;
    }
}

void isEqual_fp_test() {
    Fp::setModulo(19);

    Fp x = Fp(12);
    Fp y = Fp(12 - 19);
    Fp z = Fp(-12);

    std::cout << "[*] Fp isEq test1: ";
    if (x == y) {
        std::cout << "OK" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
    }

    std::cout << "[*] Fp isEq test2: ";
    if (x != z) {
        std::cout << "OK" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
    }

}

int main() {
    //order_test();
    //isEqual_fp_test();
    benchmark_fp();
    benchmark_ec();
}
