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
    mul(R, G, n);
    //r_mul(R, G, n);
    //montgomery_mul(R, G, n);
    //window_mul(R, G, n);

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
    mul(R, G, n);

    Point O = EC(0, 1, 0);
    if (R == O) {
        std::cout << "[*] order test: OK" << std::endl;
    } else {
        std::cout << "[*] order test: FAILED" << std::endl;
    }
}

void test_ec_muls() { // 4つのスカラー倍の計算テスト
    mpz_class p = mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    Fp::setModulo(p);

    mpz_class a = mpz_class("0", 10);
    mpz_class b = mpz_class("7", 10);

    EllipticCurve EC = EllipticCurve(a, b);
    mpz_class n = mpz_class("1154235709850086879078528375642790749043826051631415186137", 10); 
    mpz_class gx = mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_class gy = mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Point G = EC(gx, gy);
    Point R1, R2, R3, R4;
    mul(R1, G, n);
    r_mul(R2, G, n);
    montgomery_mul(R3, G, n);
    window_mul(R4, G, n);

    if (isEqual(R1,R2) && isEqual(R2,R3) && isEqual(R3,R4)) {
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

void test_GLVsecp256k1_baseMul() {
    GLV::initForsecp256k1();
    mpz_class k;
    Point R1, R2;
    k = mpz_class("3321038201388210320131380183201838214891840184028302814104802918301", 16);
    GLV::mulBase(R1, k);
    mul(R2, GLV::base, k);
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
    k = mpz_class("3321038201388210320131380183201838214891840184028302814104802918301", 16);
    mul(R, GLV::base, 382108383); 
    GLV::scalarMul(R1, R, k);
    mul(R2, R, k);
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
    Point R1, R2;
    k1 = mpz_class("11117289373161954235709850086879078528375642790749043841647", 16);
    k2 = mpz_class("DEADBEEF3921391232134374927392173937137213797392713292193", 16);
    mul(R1, GLV::base, k1 + k2); // R1 = [k+k]Base
    multipleMul(R2, GLV::base, k1, GLV::base, k2); // [k]Base + [k]Base

    std::cout << "[*] Multiple Scalar Mul test: ";
    if (R1 == R2) {
        puts("OK");
    } else {
        puts("Failed");
    }
}


int main() {
    test_GLVsecp256k1_baseMul();
    test_GLVsecp256k1_ScalarMul();
    test_MultipleScalarMul();
    test_fp_squareRoot();
    test_ECorder();
    test_ec_mul();
    test_ec_muls();
    test_isEqual_fp();

}
