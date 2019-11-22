// $ g++ mr_test.cpp -lgmpxx -lgmp -O3 -Wall -Wextra

#include<gmpxx.h>
#include<iostream>

#include<time.h>

mpz_class n, g;
mpz_class n_, r, r_; // for MR
mpz_class r2;

mpz_class X, Y, XY, S;
mpz_class T;


void setup() {
    size_t n_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    mpz_ui_pow_ui(r.get_mpz_t(), 2, n_bits); // r := n^{k+1} ( > n)
    mpz_gcdext(g.get_mpz_t(), n_.get_mpz_t(), r_.get_mpz_t(), n.get_mpz_t(), r.get_mpz_t()); // -n_*n + r_*r = 1
    mpz_mod(n_.get_mpz_t(), n_.get_mpz_t(), r.get_mpz_t()); // -n_ * n = 1 mod r
    mpz_mod(r_.get_mpz_t(), r_.get_mpz_t(), n.get_mpz_t()); // r_ * r  = 1 mod n
    if (n_ > 0) {
        mpz_sub(n_.get_mpz_t(), r.get_mpz_t(), n_.get_mpz_t()); // n_ <- -1 * n_ = r - n_
    } else {
        mpz_mul_ui(n_.get_mpz_t(), n_.get_mpz_t(), -1); // n_ <- -1 * -n_
    }
    mpz_sub_ui(T.get_mpz_t(), r.get_mpz_t(), 1); // R - 1 = 2^{k} - 1
    mpz_powm_ui(r2.get_mpz_t(), r.get_mpz_t(), 2, n.get_mpz_t()); // r2 <- r^{2} mod n
}

void MontRe(mpz_class &S, const mpz_class &XY) { // r <- XY*R^{-1}
    // S <- (XY + n * (n_ * XY mod R)) / R
    mpz_mul(S.get_mpz_t(), n_.get_mpz_t(), XY.get_mpz_t()); // n_ * XY
#if 0
    mpz_mod(S.get_mpz_t(), S.get_mpz_t(), r.get_mpz_t()); // n_ * XY mod R
#else
    mpz_and(S.get_mpz_t(), S.get_mpz_t(), T.get_mpz_t()); // mod 2^{k} = and (2^k - 1)
#endif
    mpz_mul(S.get_mpz_t(), S.get_mpz_t(), n.get_mpz_t()); // n * (n_ * XY mod R)
    mpz_add(S.get_mpz_t(), S.get_mpz_t(), XY.get_mpz_t()); // (XY + n * (n_ * XY mod R))
    mpz_tdiv_q_2exp(S.get_mpz_t(), S.get_mpz_t(), mpz_sizeinbase(n.get_mpz_t(), 2)); // (XY + n * (n_ * XY mod R)) / R

    if (S > n) {
        mpz_sub(S.get_mpz_t(), S.get_mpz_t(), n.get_mpz_t());
    }
}

void Mod(mpz_class &S, const mpz_class &XY) { // r <- XY mod n
    mpz_mod(S.get_mpz_t(), XY.get_mpz_t(), n.get_mpz_t());
}


void MMM(mpz_class &ret, const mpz_class &x, const mpz_class &y) { // r <- x*y mod n
#if 1
    mpz_class s;
    mpz_mul(s.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
    MontRe(ret, s);
    mpz_mul(s.get_mpz_t(), ret.get_mpz_t(), r2.get_mpz_t());
    MontRe(ret, s); 
#else
    mpz_class s, t;
    mpz_mul(ret.get_mpz_t(), x.get_mpz_t(), r2.get_mpz_t());
    MontRe(s, ret); // MR(x * R2)
    mpz_mul(ret.get_mpz_t(), y.get_mpz_t(), r2.get_mpz_t());
    MontRe(t, ret); // MR(y * R2)
    mpz_mul(ret.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t()); // s * t
    MontRe(s, ret); // MR(s * t)
    MontRe(ret, s); // MR(s * t)
#endif
}

int main() {
    n = mpz_class("50d45e2daa6a42f0eea586fc91fb038446fb2f2846efb70c3c308b2cf4dcb482db23c8ffe912de08d80d13e8c23a90269e23ce42a86ab8c8d6fcfe955adfc1424172928b73f495af5a7b109192564ad974b045a567b6d3fa10b7d694858af5aebf84481dffdb651c6e0f91fd96f29241a074f5819d3746c5afd1af5", 16); // n = prime(512) * prime(512)
    X = mpz_class("f5a7b109192564ad974045a567b6d3fa10b7d694858af5aebf844eea586fc91fb038446fbb10919256af5a75642f2846efb70", 16);
    Y = mpz_class("46c5a73f495af5f2846e51245adfc142417c1424172928b73f495af5a7564ad974b045a56c8ffe912de08a567b6d3fa10b7d694", 16);
    setup();

    XY = X*Y;
    Mod(S, XY);
    std::cout << S << std::endl;

    MontRe(S, XY);
    S = (S * r) % n;
    std::cout << S << std::endl;

    MMM(S, X, Y);
    std::cout << S << std::endl;


    // X*Y mod n の値が出るまでの1サイクルを計測
    const int N = 1000000;

    // Montgomery Reduction (XY * R^{-1} mod n)
    time_t begin = clock();
    for(int i = 0; i < N; i++) {
        mpz_mul(XY.get_mpz_t(), X.get_mpz_t(), Y.get_mpz_t());
        MontRe(S, XY);
        
#if 1
        mpz_mul(S.get_mpz_t(), S.get_mpz_t(), r.get_mpz_t());
        mpz_mod(S.get_mpz_t(), S.get_mpz_t(), n.get_mpz_t());
#endif
    }
    time_t end = clock();
    puts("[+] MR(X*Y)*R");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);

    // 純粋なX*Y mod n
    begin = clock();
    for(int i = 0; i < N; i++) {
        mpz_mul(XY.get_mpz_t(), X.get_mpz_t(), Y.get_mpz_t());
        mpz_mod(S.get_mpz_t(), XY.get_mpz_t(), n.get_mpz_t());
    }
    end = clock();
    puts("[+] XY \% n");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);

    // X*Y mod n の時間かかるやつ
    begin = clock();
    for(int i = 0; i < N; i++) {
        S = X*Y % n;
    }
    end = clock();
    puts("[+] XY \% n (BAD)");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);


    // MMM(X, Y) = X*Y mod n
    begin = clock();
    for(int i = 0; i < N; i++) {
        MMM(S, X, Y);
    }
    end = clock();
    puts("[+] MMM(X, Y) = MR(MR(X * Y) * R^{2})");
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / N * 1e6);


}
