#if 1
/*
    g++ -O3 -DNDEBUG -I <mimalloc>/include <mimalloc>/build/libmimalloc.a -lpthread -lgmpxx -lgmp
*/
#include <gmp.h>
#include <mimalloc.h>
static int malloc_count;

static struct UseMiMalloc {
    static void* mi_malloc_wrapper(size_t n)
    {
        //malloc_count++;
        return mi_malloc(n);
    }
     static void* mi_realloc_wrapper(void *p, size_t, size_t n)
    {
        return mi_realloc(p, n);
    }
    static void mi_free_wrapper(void *p, size_t)
    {
        mi_free(p);
    }
    UseMiMalloc()
    {
        puts("set GMP memory functions before using mpz_class");
        mp_set_memory_functions(mi_malloc_wrapper, mi_realloc_wrapper, mi_free_wrapper);
    }
    ~UseMiMalloc() {
        printf("malloc_count=%d\n", malloc_count);
    }
} g_UseMiMalloc;
#endif

#include<iostream>
#include "curve.h"
#include<time.h>

#include "FP.h"


void benchmark_ec();
void benchmark_fp();
void isEqual_fp_test();
void order_test();
void ec_mul_test();

mpz_class Fp::modulus;
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

void ec_mul_test() {
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class n = mpz_class("5792089237316195423570985008687907852837564279074904382605163141518161494337", 10); 
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC(gx, gy);
    Point R;
    R = G*n;

    Fp x, y;
    R.xy(x, y);
    mpz_class Rx = mpz_class("24468494029366207626986019034967613638108911936555812085751778627749375846788", 10);
    mpz_class Ry = mpz_class("64452411616332977820528608943388105946346351335284667932071435835634206415910", 10);
    if (x.value == Rx && y.value == Ry) {
        std::cout << "[*] EC mul test: OK" << std::endl;
    } else {
        std::cout << "[*] EC mul test: FAILED" << std::endl;
    }
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

    Point G = EC(gx, gy);
    Point R;
    R = G*n;

    //Point O = EllipticCurve::point(0, 1, 0);
    Point O = EC(0, 1, 0);
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
    order_test();
    ec_mul_test();
    isEqual_fp_test();
    benchmark_fp();
    benchmark_ec();
}
