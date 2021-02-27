#include<iostream>
#include "curve.h"
#include<time.h>

#include "FP.h"

void benchmark_scalarMul() {
    const int n = 1000;
    time_t begin, end;

    TwistedEdwardCurve::initEd();
    Point P(mpz_class("216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51a", 16), 
            mpz_class("6666666666666666666666666666666666666666666666666666666666666658", 16), 
            1);

    mpz_class x("c11f535ca9a31c6ef3ec441f75e362c1b67ae4a37121af9cd35ce110754363cd", 16);
    Point R;

    begin = clock();
    for(int i = 0; i < n; i++) {
        TwistedEdwardCurve::scalarMul(R, P, x);
    }
    end = clock();
    printf("\tEd25519 scalarMult\t time = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

int main() {
    benchmark_scalarMul();
}
