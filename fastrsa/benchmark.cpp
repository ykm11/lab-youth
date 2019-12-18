#include <gmpxx.h>
#include "rsa.h"

#include<time.h>
#include<iostream>

void benchmark_rsa() {
    mpz_class p, q, e;
    p = mpz_class("19fbd41d69aa3d86009a967db3379c63cd501f24f7", 16);
    q = mpz_class("1b6f141f98eeb619bc0360220160a5f75ea07cdf1d", 16);
    e = 65537;

    RSA rsa = RSA(p, q, e);

    mpz_class plaintext, m, c;
    plaintext = mpz_class("12345678901234567890", 10);
    rsa.encrypt(c, plaintext);

    const int n = 50000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        rsa.decrypt(m, c);
    }
    time_t end = clock();
    std::cout << "\tmod n:";
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void test_rsa() {
    mpz_class p, q, e;
    p = mpz_class("19fbd41d69aa3d86009a967db3379c63cd501f24f7", 16);
    q = mpz_class("1b6f141f98eeb619bc0360220160a5f75ea07cdf1d", 16);
    e = 65537;

    RSA rsa = RSA(p, q, e);

    mpz_class plaintext, m, c;
    plaintext = mpz_class("1844674407370955161418446744073709551614", 16);
    rsa.encrypt(c, plaintext);
    rsa.decrypt(m, c);

    std::cout << "\tmod n: ";
    if (m == plaintext) {
        puts("OK");
    } else {
        puts("FAILED");
    }
}

void benchmark_crt_rsa() {
    mpz_class p, q, e;
    p = mpz_class("19fbd41d69aa3d86009a967db3379c63cd501f24f7", 16);
    q = mpz_class("1b6f141f98eeb619bc0360220160a5f75ea07cdf1d", 16);
    e = 65537;

    RSA rsa = RSA(p, q, e);

    mpz_class plaintext, m, c;
    plaintext = mpz_class("12345678901234567890", 10);
    rsa.encrypt(c, plaintext);

    const int n = 50000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        rsa.crt_decrypt(m, c);
    }
    time_t end = clock();
    std::cout << "\tCRT:";
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void test_crt_rsa() {
    mpz_class p, q, e;
    p = mpz_class("19fbd41d69aa3d86009a967db3379c63cd501f24f7", 16);
    q = mpz_class("1b6f141f98eeb619bc0360220160a5f75ea07cdf1d", 16);
    e = 65537;

    RSA rsa = RSA(p, q, e);

    mpz_class plaintext, m, c;
    plaintext = mpz_class("1844674407370955161418446744073709551614", 16);
    rsa.encrypt(c, plaintext);
    rsa.crt_decrypt(m, c);

    std::cout << "\tCRT: ";
    if (m == plaintext) {
        puts("OK");
    } else {
        puts("FAILED");
    }
}


int main() {
    std::cout << "[*] RSA test Dec(Enc(m)) = m ?\n";
    test_rsa();
    test_crt_rsa();
    std::cout << "\n[*] Benchmark decrypt\n";
    benchmark_rsa();
    benchmark_crt_rsa();
}
