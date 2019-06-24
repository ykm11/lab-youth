#include<stdio.h>
#include<gmpxx.h>
#include<iostream>

#include "FP.h"


class Point {
public:
    mpz_class x;
    mpz_class y;
    mpz_class z;


    Point operator+(Point Q) {
        //std::cout << x << " " << y << " " << z << std::endl;
        //std::cout << Q.x << " " << Q.y << " " << Q.z << std::endl;
        Point R;
        mpz_class s, t, u, v, w, v2, v3;
        mpz_class x3, y3, z3;

        u = (Q.y*z - y*Q.z);
        v = (Q.x*z - x*Q.z);
        //u = (Q.y*z - y*Q.z) % modulus;
        //v = (Q.x*z - x*Q.z) % modulus;

        v2 = v*v;
        v3 = v2*v;
        w = u*u*z*Q.z - v3 - 2*v2*x*Q.z;
        //std::cout << v << " " << v2 << " " << v3 << " " << w << std::endl;

        x3 = v*w;
        y3 = u * (v2 * x * Q.z - w) - v3 * y * Q.z;
        z3 = v3 * z * Q.z;

        if (z3 == 0) {
            R.x = 0;
            R.y = 1;
            R.z = 0;
        } else {
            R.x = x3;
            R.y = y3;
            R.z = z3;
        }
        return R;
    }
};


class EllipticCurve {
public:
    mpz_class A;
    mpz_class B;
    mpz_class modulus;

    EllipticCurve(mpz_class a, mpz_class b, mpz_class p) {
        A = a;
        B = b;
        modulus = p;
    }
    //~EllipticCurve();


    Point operator()(mpz_class x, mpz_class y) {
        Point P;
        if ((y*y) % modulus == ((x*x*x % modulus) + (A*x % modulus) + B) % modulus) {
            P.x = x;
            P.y = y;
            P.z = 1;
        } else { // 存在しない点にはエラーを返した方が良さそう?
            P.x = 0;
            P.y = 1;
            P.z = 0;
        }
        return P; 
    }

    void printCurve() {
        printf("EllipticCurve: y^2 = x^3");
        if (A != 0) {
            printf(" + ");
            if (A < 0) {
                std::cout << (modulus + A) % modulus;
            } else {
                std::cout << A;
            }
            printf("x");
        }
        if (B != 0) {
            printf(" + ");
            if (B < 0) {
                std::cout << (modulus + B) % modulus;
            } else {
                std::cout << B;
            }
        }
        std::cout << " mod " << modulus << std::endl;
    }
};

