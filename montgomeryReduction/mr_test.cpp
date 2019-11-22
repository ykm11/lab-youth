// $ g++ mr_test.cpp -lgmpxx -lgmp -O3 -Wall -Wextra

#include "mr.h"

#include<gmpxx.h>
#include<iostream>

#include<time.h>


void test_powm() {
    mpz_class c, m, e;
    m = 1223191;
    e = 65537;

    powm(c, m, e, n);
    std::cout << "n: " << n << std::endl;
    std::cout << "pow(m, r, n): " << c << std::endl;

}

void benchmark_powm() {
    mpz_class c, m, e;
    m = 21;
    e = 65537;

    const int N = 500000;
    time_t begin = clock();
    for(int i = 0; i < N; i++) {
        powm(c, m, e, n);
    }
    time_t end = clock();
    puts("[+] pow(m, e, n)");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);
}

void benchmark_MR() {
    const int N = 1000000;
    size_t begin = clock();
    for(int i = 0; i < N; i++) {
        mpz_mul(XY.get_mpz_t(), X.get_mpz_t(), Y.get_mpz_t());
        MontRe(S, XY);
        
#if 1
        mpz_mul(S.get_mpz_t(), S.get_mpz_t(), r.get_mpz_t());
        mpz_mod(S.get_mpz_t(), S.get_mpz_t(), n.get_mpz_t());
#endif
    }
    size_t end = clock();
    puts("[+] MR(X*Y)*R");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);
}

void benchmark_MMM() {
    const int N = 1000000;
    size_t begin = clock();
    for(int i = 0; i < N; i++) {
        MMM(S, X, Y);
    }
    size_t end = clock();
    puts("[+] MMM(X, Y) = MR(MR(X * Y) * R^{2})");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);


}

void benchmark_Mod1() {
    // 純粋なX*Y mod n
    const int N = 1000000;
    size_t begin = clock();
    for(int i = 0; i < N; i++) {
        mpz_mul(XY.get_mpz_t(), X.get_mpz_t(), Y.get_mpz_t());
        mpz_mod(S.get_mpz_t(), XY.get_mpz_t(), n.get_mpz_t());
    }
    size_t end = clock();
    puts("[+] XY \% n");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);
}

void benchmark_Mod2() {
    // X*Y mod n の時間かかるやつ
    const int N = 1000000;
    size_t begin = clock();
    for(int i = 0; i < N; i++) {
        S = X*Y % n;
    }
    size_t end = clock();
    puts("[+] XY \% n (BAD)");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);
}

int main() {
    n = mpz_class("50d45e2daa6a42f0eea586fc91fb038446fb2f2846efb70c3c308b2cf4dcb482db23c8ffe912de08d80d13e8c23a90269e23ce42a86ab8c8d6fcfe955adfc1424172928b73f495af5a7b109192564ad974b045a567b6d3fa10b7d694858af5aebf84481dffdb651c6e0f91fd96f29241a074f5819d3746c5afd1af5", 16); // n = prime(512) * prime(512)
    X = mpz_class("f5a7b109192564ad974045a567b6d3fa10b7d694858af5aebf844eea586fc91fb038446fbb10919256af5a75642f2846efb70", 16);
    Y = mpz_class("46c5a73f495af5f2846e51245adfc142417c1424172928b73f495af5a7564ad974b045a56c8ffe912de08a567b6d3fa10b7d694", 16);
    setup();

    XY = X*Y;

    //test_powm();
    benchmark_powm();
    return 0;
    
    // X*Y mod n の値が出るまでの1サイクルを計測

    // Montgomery Reduction (XY * R^{-1} mod n)
    benchmark_MR();
    
    // MMM(X, Y) = X*Y mod n
    benchmark_MMM();

    benchmark_Mod1();
    benchmark_Mod2();

}
