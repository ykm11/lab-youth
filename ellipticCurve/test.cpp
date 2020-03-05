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


void test_ec_mul() {
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
    naf_mul(R, G, n);

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

void test_ec_sub() {
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class n = mpz_class("5792089237316195423570985008687907852837564279074904382605163141518161494337", 10); 
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC(gx, gy);
    Point R5, R3, R2;
    r_mul(R2, G, 2);
    r_mul(R3, G, 3);
    r_mul(R5, G, 5);

    std::cout << "[*] EC sub test: ";
    if ((R5 - R2) == R3) {
        std::cout << "OK\n";
    } else {
        std::cout << "FAILED\n";
    }
}

void test_jacobi_ec_dbl() {
    /*
    p = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
    a = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC
    b = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B

    Gx = 6B17D1F2 E12C4247 F8BCE6E5 63A440F2 77037D81 2DEB33A0 F4A13945 D898C296 
    Gy = 4FE342E2 FE1A7F9B 8EE7EB4A 7C0F9E16 2BCE3357 6B315ECE CBB64068 37BF51F5
    */
    mpz_class p = mpz_class("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16);
    Fp::setModulo(p);

    mpz_class a = mpz_class("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC", 16);
    mpz_class b = mpz_class("5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B", 16);
    EllipticCurve EC = EllipticCurve(a, b);

    mpz_class gx = mpz_class("6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
    mpz_class gy = mpz_class("4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);
    jPoint G, G2;
    G = jPoint(gx, gy, 1);
    EllipticCurve::dbl(G2, G);

    Point P, P2;
    P = Point(gx, gy, 1);
    EllipticCurve::dbl(P2, P);

    Fp u, v, s, t;
    G2.xy(u, v);
    P2.xy(s, t);
    std::cout << "[*] EC(jacobi) dbl test (a = -3): ";
    if (u == s && v == t) {
        std::cout << "OK\n";
    } else {
        std::cout << "FAILED\n";
    }

    p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    a = mpz_class("0", 10);
    b = mpz_class("7", 10);

    EC = EllipticCurve(a, b);
    gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    G = jPoint(gx, gy, 1);
    EllipticCurve::dbl(G2, G);

    P = Point(gx, gy, 1);
    EllipticCurve::dbl(P2, P);

    G2.xy(u, v);
    P2.xy(s, t);
    std::cout << "[*] EC(jacobi) dbl test (a != -3): ";
    if (u == s && v == t) {
        std::cout << "OK\n";
    } else {
        std::cout << "FAILED\n";
    }
}

void test_jacobi_ec_add() {
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class n = mpz_class("5792089237316195423570985008687907852837564279074904382605163141518161494337", 10); 
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    jPoint G1 = jPoint(gx, gy, 1);

    Point P = EC(gx, gy);
    Point P2;
    EllipticCurve::dbl(P2, P);
    Fp px, py;
    P2.xy(px, py);
    jPoint G2;
    EllipticCurve::dbl(G2, G1);

    jPoint G3, G4;
    add(G3, G1, G2);
    EllipticCurve::dbl(G4, G2);

    Point P3, P4;
    add(P3, P2, P);
    add(P4, P3, P);

    Fp u, v, s, t;
    G4.xy(u, v);
    P4.xy(s, t);

    std::cout << "[*] EC(jacobi) add test: ";
    if (u == s && v == t) {
        std::cout << "OK\n";
    } else {
        std::cout << "FAILED\n";
    }
}

void test_ECorder() {
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
    r_mul(R, G, n);

    Point O = EC(0, 1, 0);
    if (R == O) {
        std::cout << "[*] order test: OK" << std::endl;
    } else {
        std::cout << "[*] order test: FAILED" << std::endl;
    }
}

void test_ec_muls() { 
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class n = mpz_class("1154235709850086879078528375642790749043826051631415186137", 10); 
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC(gx, gy);
    Point R1, R2, R3, R4, R5;
    l_mul(R1, G, n);
    r_mul(R2, G, n);
    montgomery_mul(R3, G, n);
    window_mul(R4, G, n);
    naf_mul(R5, G, n);

    if (isEqual(R1,R2) && isEqual(R2,R3) && isEqual(R3,R4) && isEqual(R4, R5)) {
        std::cout << "[*] EC muls test: OK" << std::endl;
    } else {
        std::cout << "[*] EC muls test: FAILED" << std::endl;
    }
}

void test_isEqual_fp() {
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

void test_fp_squareRoot() {
    std::cout << "[*] Fp squareRoot test: ";
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);
    Fp x, r;

    r.value = 0;
    for (int i = 0; i < 100; i++) {
        mpz_random(x.value.get_mpz_t(), 4);
        if (Fp::squareRoot(r, x)) {
            if (x.value != (r.value*r.value % Fp::modulus)) {
                std::cout << "Failed\n";
                return;
            }
        }
    }
    std::cout << "OK\n";
}

void test_GLV_decomposing() {
    std::cout << "[*] GLV decomposing-k secp256k1 test: ";
    GLV::initForsecp256k1();
    mpz_class k = mpz_class("68db8bac710cb295e9e1b089a0275253db6a70997889e1c902cb5018e8bd5", 16);
    mpz_class k0, k1;

    GLV::decomposing_k(k0, k1, k);
    if ((k0 + k1*GLV::lmd) % GLV::order == k % GLV::order) {
        std::cout << "OK\n";
    } else {
        std::cout << "Failed\n";
    }
}

void test_GLVsecp256k1_baseMul() {
    GLV::initForsecp256k1();
    mpz_class k;
    Point R1, R2;
    k = mpz_class("68db8bac710cb295e9e1b089a0275253db6a70997889e1c902cb5018e8bd5", 16);
    GLV::mulBase(R1, k);
    r_mul(R2, GLV::base, k);
    std::cout << "[*] GLV base mul test: ";
    if (R1 == R2) {
        puts("OK");
    } else {
        puts("Failed");
    }
}

void test_GLVsecp256k1_ScalarMul() {
    GLV::initForsecp256k1();
    mpz_class k;
    Point R, R1, R2;
    k = mpz_class("68db8bac710cb295e9e1b089a0275253db6a70997889e1c902cb5018e8bd5", 16);
    r_mul(R, GLV::base, 382108383); 
    GLV::scalarMul(R1, R, k);
    r_mul(R2, R, k);
    std::cout << "[*] GLV scalar mul test: ";
    if (R1 == R2) {
        puts("OK");
    } else {
        puts("Failed");
    }
}

void test_MultipleScalarMul() {
    GLV::initForsecp256k1();
    mpz_class k1, k2;
    Point R, R1, R2;
    k1 = mpz_class("11117289373161954235709850086879078528375642790749043841647", 16);
    k2 = mpz_class("DEADBEEF3921391232134374927392173937137213797392713292193", 16);
    GLV::mulBase(R, 38210831383); 
    r_mul(R1, R, k1 + k2); // R1 = [k+k]Base
    multipleMul(R2, R, k1, R, k2); // [k]Base + [k]Base

    std::cout << "[*] Multiple Scalar Mul test: ";
    if (R1 == R2) {
        puts("OK");
    } else {
        puts("Failed");
    }
}


int main() {
    test_GLV_decomposing();
    test_GLVsecp256k1_baseMul();
    test_GLVsecp256k1_ScalarMul();
    test_MultipleScalarMul();
    test_fp_squareRoot();
    test_ECorder();
    test_ec_sub();
    test_ec_mul();
    test_ec_muls();
    test_isEqual_fp();

    test_jacobi_ec_add();
    test_jacobi_ec_dbl();

}
