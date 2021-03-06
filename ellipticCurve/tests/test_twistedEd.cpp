#include "curve.h"
#include "FP.h"
#include <iostream>
#include <string.h>

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

void test_TwistedEd_baseMult() {
    TwistedEdwardCurve::initEd();
    mpz_class x("c11f535ca9a31c6ef3ec441f75e362c1b67ae4a37121af9cd35ce110754363cd", 16);
    Point R;
    TwistedEdwardCurve::baseMult(R, x);

    Point correct_R(mpz_class("7e2eecc8c6ef10e861af05ac15e1361f4f7e15d2023a1617f5f0cb10cea940e8", 16), 
            mpz_class("6c8be36d9349bfaae4d09dfb379b49417dc65afbf41e9815a14be0d53abdfe88", 16), 
            1);

    printf("[+] Ed25519 baseMult TEST: ");
    if (R == correct_R) {
        puts("OK");
    } else {
        puts("FAILED");
    }
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

void test_Ed25519_encodePoint() {
    uint8_t buf[32];
    TwistedEdwardCurve::initEd();
    Point R;

    mpz_class x("43145378195037189031758314739173849104", 16);
    TwistedEdwardCurve::scalarMul(R, TwistedEdwardCurve::base_, x);
    TwistedEdwardCurve::encodePoint(buf, R);

    uint8_t act_buf[] = {
        0xd7, 0xa7, 0x24, 0xa4, 0x00, 0xd2, 0x46, 0x4a, 
        0x8a, 0xc3, 0x2d, 0x18, 0x7d, 0x30, 0x53, 0xe4, 
        0x62, 0x58, 0x04, 0x86, 0x86, 0x92, 0xaa, 0x44, 
        0x47, 0x76, 0x9c, 0x52, 0x7c, 0xc0, 0x4f, 0x32,
    };

    printf("[+] encodePoint TEST: ");
    if (memcmp(buf, act_buf, 32) == 0) {
        puts("OK");
    } else {
        puts("FAILED");
    }
}

void test_Ed25519_xrecover() {
    Fp rx, ry;
    Fp x;
    Point R;
    TwistedEdwardCurve::initEd();

    R = TwistedEdwardCurve::base_;
    R.xy(rx, ry);
    TwistedEdwardCurve::xrecover(x, ry);

    dump(x);
    dump(rx);
    TwistedEdwardCurve::baseMult(R, 100);
    R.xy(rx, ry);
    TwistedEdwardCurve::xrecover(x, ry);
}


int main() {
    //test_TwistedEd_add();
    test_Ed25519_order();
    test_TwistedEd_scalarMul();
    test_TwistedEd_baseMult();
    test_Ed25519_encodePoint();
    test_Ed25519_xrecover();
}
