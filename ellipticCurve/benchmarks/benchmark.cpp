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

void benchmark_ec_add_secp256k1() {
    std::cout << "[*] secp256k1 add benchmark\n";
    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC.point(gx, gy);
    Point R, Q;
    EllipticCurve::dbl(R, G);

    const int n = 100000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        add(Q, R, G);
    }
    time_t end = clock();
    printf("\tProjection\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    jPoint G1(gx, gy, 1);
    jPoint R1, Q1;
    EllipticCurve::dbl(R1, G1);

    begin = clock();
    for(int i = 0; i < n; i++) {
        add(Q1, R1, G1);
    }
    end = clock();
    printf("\tJacobian\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_ec_add_P256() {
    std::cout << "[*] NIST P256 add benchmark\n";
    mpz_class a = mpz_class("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC", 16);
    mpz_class b = mpz_class("5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B", 16);
    mpz_class p = mpz_class("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16);
    Fp::setModulo(p);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class gx = mpz_class("6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
    mpz_class gy = mpz_class("4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);

    Point G = EC.point(gx, gy);
    Point R, Q;
    EllipticCurve::dbl(R, G);

    const int n = 100000;
    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        add(Q, R, G);
    }
    time_t end = clock();
    printf("\tProjection\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    jPoint G1(gx, gy, 1);
    jPoint R1, Q1;
    EllipticCurve::dbl(R1, G1);

    begin = clock();
    for(int i = 0; i < n; i++) {
        add(Q1, R1, G1);
    }
    end = clock();
    printf("\tJacobian\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_ec_dbl_secp256k1() {
    std::cout << "[*] secp256k1 dbl benchmark\n";
    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = Point(gx, gy, 1);
    Point R = Point(0, 1, 0);
    EllipticCurve::dbl(R, G);
    const int n = 100000;

    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        EllipticCurve::dbl(R, R);
    }
    time_t end = clock();
    printf("\tProjection\t time = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    jPoint jG = jPoint(gx, gy, 1);
    jPoint jR;
    EllipticCurve::dbl(jR, jG);
    begin = clock();
    for(int i = 0; i < n; i++) {
        EllipticCurve::dbl(jR, jR);
    }
    end = clock();
    printf("\tJacobian\t time = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_ec_dbl_P256() {
    std::cout << "[*] NIST P256 dbl benchmark\n";
    mpz_class p = mpz_class("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16);
    mpz_class a = mpz_class("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC", 16);
    mpz_class b = mpz_class("5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B", 16);
    Fp::setModulo(p);
    EllipticCurve EC = EllipticCurve(a, b);

    mpz_class gx = mpz_class("6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
    mpz_class gy = mpz_class("4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);

    Point G = Point(gx, gy, 1);
    Point R = Point(0, 1, 0);
    EllipticCurve::dbl(R, G);
    const int n = 100000;

    time_t begin = clock();
    for(int i = 0; i < n; i++) {
        EllipticCurve::dbl(R, R);
    }
    time_t end = clock();
    printf("\tProjection\t time = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    jPoint jG = jPoint(gx, gy, 1);
    jPoint jR;
    EllipticCurve::dbl(jR, jG);
    begin = clock();
    for(int i = 0; i < n; i++) {
        EllipticCurve::dbl(jR, jR);
    }
    end = clock();
    printf("\tJacobian\t time = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}


void benchmark_ec_mul_secp256k1() {
    std::cout << "[*] secp256k1 mul benchmark\n";
    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class q = mpz_class("35413290456945547306056027344241947654113124042893875504030834988367514449923", 10);
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC.point(gx, gy);
    Point R;
    const int n = 1000;
    time_t begin, end;

    std::cout << "\t[Proj]\n";
    // 左向きバイナリ法
    std::cout << "\tRtL Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        l_mul(R, G, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    // 右向きバイナリ法
    std::cout << "\tLtR Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        r_mul(R, G, q);
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

    std::cout << "\t[Jacobi]\n";
    jPoint G1, R1;
    G1 = jPoint(gx, gy, 1);
    std::cout << "\tRtL Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        l_mul(R1, G1, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\tLtR Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        r_mul(R1, G1, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\twin-sli(w=2)";
    begin = clock();
    for(int i = 0; i < n; i++) {
        window_mul(R1, G1, q);
    }
    end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\tNaf";
    begin = clock();
    for(int i = 0; i < n; i++) {
        naf_mul(R1, G1, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}


void benchmark_ec_mul_P256() {
    std::cout << "[*] NIST P256 mul benchmark\n";
    mpz_class p = mpz_class("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16);
    mpz_class a = mpz_class("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC", 16);
    mpz_class b = mpz_class("5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B", 16);
    Fp::setModulo(p);
    EllipticCurve EC = EllipticCurve(a, b);

    mpz_class gx = mpz_class("6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
    mpz_class gy = mpz_class("4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);

    mpz_class q = mpz_class("35413290456945547306056027344241947654113124042893875504030834988367514449923", 10);

    Point G = Point(gx, gy, 1);
    Point R;
    const int n = 1000;
    time_t begin, end;

    std::cout << "\t[Proj]\n";
    // 左向きバイナリ法
    std::cout << "\tRtL Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        l_mul(R, G, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    // 右向きバイナリ法
    std::cout << "\tLtR Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        r_mul(R, G, q);
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

    std::cout << "\t[Jacobi]\n";
    jPoint G1, R1;
    G1 = jPoint(gx, gy, 1);
    std::cout << "\tRtL Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        l_mul(R1, G1, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\tLtR Bin";
    begin = clock();
    for(int i = 0; i < n; i++) {
        r_mul(R1, G1, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\twin-sli(w=2)";
    begin = clock();
    for(int i = 0; i < n; i++) {
        window_mul(R1, G1, q);
    }
    end = clock();
    printf("\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

    std::cout << "\tNaf";
    begin = clock();
    for(int i = 0; i < n; i++) {
        naf_mul(R1, G1, q);
    }
    end = clock();
    printf("\t\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}


void benchmark_fp() {
    std::cout << "[*] Fp invmod benchmark\n";
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);
    Fp x(mpz_class("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10));
    Fp y;

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

    x = Fp(mpz_class("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10));

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

    x = Fp(mpz_class("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10)); 

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
    //benchmark_fp_sqareRoot();
    //benchmark_fp();
    //benchmark_sqr();
    benchmark_ec_add_secp256k1();
    benchmark_ec_dbl_secp256k1();
    benchmark_ec_add_P256();
    benchmark_ec_dbl_P256();
    benchmark_ec_mul_secp256k1();
    benchmark_ec_mul_P256();

    benchmark_GLV_decomposing();
    benchmark_GLVbaseMul(); 
    benchmark_MultipleScalarMul();
}
