#include "curve.h"
#include "FP.h"
#include <iostream>

#include <gmpxx.h>


#define SECP521_SIZE 9

/*

b =     0x0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00

Gx =    0x00C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66
Gy =    0x011839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650

n = 0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409

*/
void test_FP521sqr() {
    mp_limb_t p[SECP521_SIZE] = {
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };
    Fp::setModulo(p);

    mp_limb_t a[SECP521_SIZE] = {
        0x2137201477420138L, 0x6313232130472104L, 0x3639164821638291L, 0x8371372917431694L,
        0x1739278732747374L, 0x1732913789724789L, 0x7296392641693692L, 0x32dead13L,
    };

    mp_limb_t c[SECP521_SIZE] = {
        0xb54ed68553c0d66aL, 0x472128ceb6f2fe70L, 0xe1c253e06047ec35L, 0xad4178b877c07038L,
        0x45d683541d9077b8L, 0x4576d60c33560266L, 0xe15fd51f6649f17eL, 0xf24e2cd520844c75L,
        0x119L,
    };

    Fp x(a);
    Fp z;
    sqr(z, x);

    Fp act_z(c);
    std::cout << "[*] Fp sqr test: ";
    if (z == act_z) {
        std::cout << "OK\n";
    } else {
        std::cout << "not OK\n";
    }
}


void test_FP521mul() {
    mp_limb_t p[SECP521_SIZE] = {
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };
    Fp::setModulo(p);

    mp_limb_t a[SECP521_SIZE] = {
        0x2137201477420138L, 0x6313232130472104L, 0x3639164821638291L, 0x8371372917431694L,
        0x1739278732747374L, 0x1732913789724789L, 0x7296392641693692L, 0x32dead13L,
    };

    mp_limb_t b[SECP521_SIZE] = {
        0x8032810382412080L, 0x2427104721732301L, 0x3281037402174017L, 0x7201839284217487L, 
        0x1737201372138214L, 0x3747217427103782L, 0x13729dee74L,
    };

    mp_limb_t c[SECP521_SIZE] = {
        0xd0703a58517ed36aL, 0xf528ff505bb4f563L, 0x4d076ef4b10b90b8L, 0x716bfdbf3276c02fL,
        0xd352a66d881d9a7dL, 0x5de5b5a08d52fdf8L, 0xcf56240ba7dabaa3L, 0x823193d5094ac774L,
        0x183L,
    };

    Fp x(a);
    Fp y(b);
    Fp z;
    mul(z, x, y);

    Fp act_z(c);
    std::cout << "[*] Fp mul test: ";
    if (z == act_z) {
        std::cout << "OK\n";
    } else {
        std::cout << "not OK\n";
    }
}


void test_secp521r() {
    mp_limb_t p[SECP521_SIZE] = {
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };
    mp_limb_t a[SECP521_SIZE] = {
        0xfffffffffffffffc, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };

    mp_limb_t b[SECP521_SIZE] = {
        0xef451fd46b503f00L, 0x3573df883d2c34f1L, 0x1652c0bd3bb1bf07L, 0x56193951ec7e937bL,
        0xb8b489918ef109e1L, 0xa2da725b99b315f3L, 0x929a21a0b68540eeL, 0x953eb9618e1c9a1fL,
        0x51L,
    };
    mp_limb_t gx[SECP521_SIZE] = {
        0xf97e7e31c2e5bd66L, 0x3348b3c1856a429bL, 0xfe1dc127a2ffa8deL, 0xa14b5e77efe75928L,
        0xf828af606b4d3dbaL, 0x9c648139053fb521L, 0x9e3ecb662395b442L, 0x858e06b70404e9cdL,
        0xc6L,
    };
    mp_limb_t gy[SECP521_SIZE] = {
        0x88be94769fd16650L, 0x353c7086a272c240L, 0xc550b9013fad0761L, 0x97ee72995ef42640L,
        0x17afbd17273e662cL, 0x98f54449579b4468L, 0x5c8a5fb42c7d1bd9L, 0x39296a789a3bc004L,
        0x118L,
    };
    Fp::setModulo(p);

    Fp aa(a);
    Fp bb(b);
    EllipticCurve E(aa, bb);

    Fp Gx(gx);
    Fp Gy(gy);
    Fp Gz(1);
    jPoint G(Gx, Gy, Gz); 

    mpz_class q("14eba7e7ab32edc54a71f83be782d4d4be311be0d4e6c0b80bd36dffba1212020f7faacb3e839b5a2dca813576beb09d4d43db70f8b01bc1dfa18d18fde7cfd7147", 16);
    jPoint R;
    naf_mul(R, G, q);

    mpz_class rx("1c6cbe47fbd79919e77522951800c04b519c10dacd1cc1e8e35f79a1cdca44b24a2f0bab26c5348609a228820cb5842095bf9fd518296941d1c0fc734cda313fb81", 16);
    mpz_class ry("bc994a737d1c40174c56661e7481e569ddab7c7b23f621a161b6930f6006d98712337129c96bc4219e9b964878726af0108211e7414617911f3097ca89d3092bf9", 16);

    Fp Rx(rx);
    Fp Ry(ry);
    jPoint act_R(Rx, Ry, Gz);
    std::cout << "[*] secp521r1 mul test: ";
    if (isEqual(R, act_R)) {
        std::cout << "OK\n";
    } else {
        std::cout << "Failed\n";
    }
}

int main() {
    test_FP521mul();
    test_FP521sqr();
    test_secp521r();
}
