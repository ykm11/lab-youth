#ifdef USE_MIMALLOC
/*
    g++ -O3 -DNDEBUG -I <mimalloc>/include <mimalloc>/build/libmimalloc.a -lpthread -lgmpxx -lgmp
*/
#include <gmp.h>
#include <mimalloc.h>
static int malloc_count;

static struct UseMiMalloc {
    static void* mi_malloc_wrapper(size_t n)
    {
        malloc_count++;
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

void benchmark_ec_add();
void benchmark_ec_mul();
void benchmark_fp();

void benchmark_ec_add() {
    std::cout << "[*] EC add benchmark\n";
    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC.point(gx, gy);
    Point R = Point(0, 1, 0);
    const int n = 100000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        add(R, R, G);
    }
    time_t end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_ec_dbl() {
    std::cout << "[*] EC dbl benchmark\n";
    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class q = mpz_class("117289373161954235709850086879078528375642790749043841647", 10);
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC.point(gx, gy);
    Point R = Point(0, 1, 0);
    EllipticCurve::dbl(R, G);
    const int n = 100000;

    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        EllipticCurve::dbl(R, R);
    }
    time_t end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}


void benchmark_ec_mul() {
    std::cout << "[*] EC mul benchmark\n";
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
    const int n = 100;
    time_t begin, end;

    std::cout << "\tRtL Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        mul(R, G, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\tLtR Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        r_mul(R, G, q);
        //montgomery_mul(R, G, q);
        //window_mul(R, G, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\twin-sli(w=2)";
    begin = clock();
    for(int i = 0; i < n; i++) {
        window_mul(R, G, q);
    }
    end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\tNaf";
    begin = clock();
    for(int i = 0; i < n; i++) {
        naf_mul(R, G, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}


void benchmark_fp() {
    std::cout << "[*] Fp invmod benchmark\n";
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);
    Fp x, y;

    x = Fp("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10); 
    const int n = 100000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        invmod(y, x); // x <- y^{-1}
    }
    time_t end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}


void benchmark_sqr() {
    std::cout << "[*] sqr benchmark\n";
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);
    Fp x, r;

    x = Fp("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10); 

    const int n = 1000000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        sqr(r, x);
        //mul(r, x, x);
    }
    time_t end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_fp_sqareRoot() {
    std::cout << "[*] Fp squareRoot benchmark\n";
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);
    Fp x, r;

    x = Fp("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10); 

    const int n = 10000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        Fp::squareRoot(r, x);
    }
    time_t end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_GLV_decomposing() {
    std::cout << "[*] GLV decomposing-k secp256k1 benchmark\n";
    GLV::initForsecp256k1();
    mpz_class k = mpz_class("68db8bac710cb295e9e1b089a0275253db6a70997889e1c902cb5018e8bd5", 16);
    mpz_class k0, k1;
    const int n = 10000;

    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        GLV::decomposing_k(k0, k1, k);
    }
    time_t end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_GLVbaseMul() {
    std::cout << "[*] secp256k1 base mul benchmark\n";
    GLV::initForsecp256k1();
    Point R;
    mpz_class k = mpz_class("68db8bac710cb295e9e1b089a0275253db6a70997889e1c902cb5018e8bd5", 16);
    const int n = 1000;

    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        GLV::mulBase(R, k);
    }
    time_t end = clock();
    std::cout << "\tGLV";
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    begin = clock();
    for(int i = 0; i < n; i++) {
        naf_mul(R, GLV::base, k);
    }
    end = clock();
    std::cout << "\tUsual";
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_MultipleScalarMul() {
    std::cout << "[*] MultipleMul v.s Mul\n";
    GLV::initForsecp256k1();
    Point R;
    mpz_class k1 = mpz_class("11117289373161954235709850086879078528375642790749043841647", 16);
    mpz_class k2 = mpz_class("DEADBEEF3921391232134374927392173937137213797392713292193", 16);
    const int n = 1000;

    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        multipleMul(R, GLV::base, k1, GLV::base, k2);
    }
    time_t end = clock();
    std::cout << "\t[k1]Base+[k2]Base";
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    mpz_class k = k1+k2;
    begin = clock();
    for(int i = 0; i < n; i++) {
        naf_mul(R, GLV::base, k); 
    }
    end = clock();
    std::cout << "\t[k1+k2]Base";
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}


int main() {
    benchmark_fp_sqareRoot();
    benchmark_fp();
    benchmark_sqr();
    benchmark_ec_add();
    benchmark_ec_dbl();
    benchmark_ec_mul();

    benchmark_GLV_decomposing();
    benchmark_GLVbaseMul(); 
    benchmark_MultipleScalarMul();
}
