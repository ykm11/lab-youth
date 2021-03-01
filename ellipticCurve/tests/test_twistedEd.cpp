#include "curve.h"
#include "FP.h"
#include <iostream>

#include <gmpxx.h>

#define SECP521_SIZE 9

/*
http://ed25519.cr.yp.to/python/ed25519.py

d = 0x52036cee2b6ffe738cc740797779e89800700a4d4141d8ab75eb4dca135978a3L
Bx = 0x216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51aL
By = 0x6666666666666666666666666666666666666666666666666666666666666658L

order = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed

*/

void test_TwistedEd_add() {
    TwistedEdwardCurve::initEd();
    Point P(mpz_class("216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51a", 16), 
            mpz_class("6666666666666666666666666666666666666666666666666666666666666658", 16), 
            1);

    Point Q;
    TwistedEdwardCurve::Pdbl(Q, P);
    dump(Q);

    TwistedEdwardCurve::Padd(Q, P, P);
    dump(Q);
}

void test_TwistedEd_scalarMul() {
    TwistedEdwardCurve::initEd();
    Point P(mpz_class("216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51a", 16), 
            mpz_class("6666666666666666666666666666666666666666666666666666666666666658", 16), 
            1);

    mpz_class x("c11f535ca9a31c6ef3ec441f75e362c1b67ae4a37121af9cd35ce110754363cd", 16);
    Point R;
    TwistedEdwardCurve::scalarMul(R, P, x);

    Point correct_R(mpz_class("7e2eecc8c6ef10e861af05ac15e1361f4f7e15d2023a1617f5f0cb10cea940e8", 16), 
            mpz_class("6c8be36d9349bfaae4d09dfb379b49417dc65afbf41e9815a14be0d53abdfe88", 16), 
            1);

    printf("[+] Ed25519 scalarMult TEST: ");
    if (R == correct_R) {
        puts("OK");
    } else {
        puts("FAILED");
    }
}

void test_Ed25519_order() {
    mpz_class x("1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed", 16);
    TwistedEdwardCurve::initEd();
    Point P(mpz_class("216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51a", 16), 
            mpz_class("6666666666666666666666666666666666666666666666666666666666666658", 16), 
            1);
    Point R;
    TwistedEdwardCurve::scalarMul(R, P, x);
    dump(R);
}

int main() {
    //test_TwistedEd_add();
    test_Ed25519_order();
    test_TwistedEd_scalarMul();
}
