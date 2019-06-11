#include<stdio.h>
#include<gmpxx.h>
#include<iostream>


struct point {
    mpz_class x;
    mpz_class y;
    bool isUnit;
};

struct EllipticCurve {
    mpz_class A;
    mpz_class B;
    mpz_class modulus;

    void printCurve() {
        printf("EllipticCurve: y^2 = x^3 + ");
        if (A < 0) {
            std::cout << (modulus + A) % modulus;
        } else {
            std::cout << A;
        }
        printf("x + ");
        if (B < 0) {
            std::cout << (modulus + B) % modulus;
        } else {
            std::cout << B;
        }
        std::cout << " mod " << modulus << std::endl;
    }

    point Point(mpz_class x, mpz_class y) {
        point P;
        if ((y*y) % modulus == ((x*x*x % modulus) + (A*x % modulus) + B) % modulus) {
            P.x = x;
            P.y = y;
            P.isUnit = false;
        } else {
            P.x = 0;
            P.y = 1;
            P.isUnit = true;
        }
        return P;
    }
};

EllipticCurve NewCurve(mpz_class A, mpz_class B, mpz_class modulus) {
    EllipticCurve EC;
    EC.A = A % modulus;
    EC.B = B % modulus;
    EC.modulus = modulus;
    return EC; 
}

