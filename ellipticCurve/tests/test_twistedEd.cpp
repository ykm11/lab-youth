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

*/

void test_TwistedEd_add() {
    TwistedEdwardCurve::initEd();
    Point P(mpz_class("216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51a", 16), 
            mpz_class("6666666666666666666666666666666666666666666666666666666666666658", 16), 
            1);
    dump(P);

    Point Q;
    TwistedEdwardCurve::Pdbl(Q, P);
    dump(Q);

    TwistedEdwardCurve::Padd(Q, P, P);
    dump(Q);
}


int main() {
    test_TwistedEd_add();
}
