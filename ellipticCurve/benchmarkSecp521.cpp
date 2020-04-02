#ifndef SECP521
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

#include "curve.h"
#include "FP.h"
#include <iostream>

#include <time.h>

/*

b =     0x0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00

Gx =    0x00C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66
Gy =    0x011839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650

n = 0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409

*/
void benchmark_FP521sqr() {
#ifdef SECP521
    uint64_t p[SIZE] = {
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };

    uint64_t a[SIZE] = {
        0x2137201477420138L, 0x6313232130472104L, 0x3639164821638291L, 0x8371372917431694L,
        0x1739278732747374L, 0x1732913789724789L, 0x7296392641693692L, 0x32dead13L,
    };
#else
    mpz_class p("1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    mpz_class a("32DEAD137296392641693692173291378972478917392787327473748371372917431694363916482163829163132321304721042137201477420138", 16);
#endif
    Fp::setModulo(p);
    Fp x(a);
    Fp z;

    const int n = 100000;
    time_t begin, end;
    begin = clock();
    for (int i = 0; i < n; i++) {
        sqr(z, x);
    }
    end = clock();
    printf("\tFp 521 sqr \ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_FP521mul() {
#ifdef SECP521
    uint64_t p[SIZE] = {
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };

    uint64_t a[SIZE] = {
        0x2137201477420138L, 0x6313232130472104L, 0x3639164821638291L, 0x8371372917431694L,
        0x1739278732747374L, 0x1732913789724789L, 0x7296392641693692L, 0x32dead13L,
    };

    uint64_t b[SIZE] = {
        0x8032810382412080L, 0x2427104721732301L, 0x3281037402174017L, 0x7201839284217487L, 
        0x1737201372138214L, 0x3747217427103782L, 0x13729dee74L,
    };
#else
    mpz_class p("1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    mpz_class a("32DEAD137296392641693692173291378972478917392787327473748371372917431694363916482163829163132321304721042137201477420138", 16);
    mpz_class b("13729DEE74374721742710378217372013721382147201839284217487328103740217401724271047217323018032810382412080", 16);
#endif
    Fp::setModulo(p);
    Fp x(a);
    Fp y(b);
    Fp z;

    const int n = 100000;
    time_t begin, end;
    begin = clock();
    for (int i = 0; i < n; i++) {
        mul(z, x, y);
    }
    end = clock();
    printf("\tFp 521 mul \ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);

}

void benchmark_secp521r() {
#ifdef SECP521
    uint64_t p[SIZE] = {
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };
    uint64_t a[SIZE] = {
        0xfffffffffffffffc, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };

    uint64_t b[SIZE] = {
        0xef451fd46b503f00L, 0x3573df883d2c34f1L, 0x1652c0bd3bb1bf07L, 0x56193951ec7e937bL,
        0xb8b489918ef109e1L, 0xa2da725b99b315f3L, 0x929a21a0b68540eeL, 0x953eb9618e1c9a1fL,
        0x51L,
    };
    uint64_t gx[SIZE] = {
        0xf97e7e31c2e5bd66L, 0x3348b3c1856a429bL, 0xfe1dc127a2ffa8deL, 0xa14b5e77efe75928L,
        0xf828af606b4d3dbaL, 0x9c648139053fb521L, 0x9e3ecb662395b442L, 0x858e06b70404e9cdL,
        0xc6L,
    };
    uint64_t gy[SIZE] = {
        0x88be94769fd16650L, 0x353c7086a272c240L, 0xc550b9013fad0761L, 0x97ee72995ef42640L,
        0x17afbd17273e662cL, 0x98f54449579b4468L, 0x5c8a5fb42c7d1bd9L, 0x39296a789a3bc004L,
        0x118L,
    };
#else
    mpz_class p("1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    mpz_class a("1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffc", 16);
    mpz_class b("0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00", 16);
    mpz_class gx("00C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66", 16);
    mpz_class gy("011839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650", 16);
#endif

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

    const int n = 1000;
    time_t begin, end;

    begin = clock();
    for (int i = 0; i < n; i++) {
        naf_mul(R, G, q);
    }
    end = clock();
    printf("\tJacobi mul\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark_secp521r_dbl() {
#ifdef SECP521
    uint64_t p[SIZE] = {
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };
    uint64_t a[SIZE] = {
        0xfffffffffffffffc, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };

    uint64_t b[SIZE] = {
        0xef451fd46b503f00L, 0x3573df883d2c34f1L, 0x1652c0bd3bb1bf07L, 0x56193951ec7e937bL,
        0xb8b489918ef109e1L, 0xa2da725b99b315f3L, 0x929a21a0b68540eeL, 0x953eb9618e1c9a1fL,
        0x51L,
    };
    uint64_t gx[SIZE] = {
        0xf97e7e31c2e5bd66L, 0x3348b3c1856a429bL, 0xfe1dc127a2ffa8deL, 0xa14b5e77efe75928L,
        0xf828af606b4d3dbaL, 0x9c648139053fb521L, 0x9e3ecb662395b442L, 0x858e06b70404e9cdL,
        0xc6L,
    };
    uint64_t gy[SIZE] = {
        0x88be94769fd16650L, 0x353c7086a272c240L, 0xc550b9013fad0761L, 0x97ee72995ef42640L,
        0x17afbd17273e662cL, 0x98f54449579b4468L, 0x5c8a5fb42c7d1bd9L, 0x39296a789a3bc004L,
        0x118L,
    };
#else
    mpz_class p("1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    mpz_class a("1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffc", 16);
    mpz_class b("0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00", 16);
    mpz_class gx("00C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66", 16);
    mpz_class gy("011839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650", 16);
#endif

    Fp::setModulo(p);

    Fp aa(a);
    Fp bb(b);
    EllipticCurve E(aa, bb);

    Fp Gx(gx);
    Fp Gy(gy);
    Fp Gz(1);
    jPoint G(Gx, Gy, Gz); 
    jPoint R;

    const int n = 100000;
    time_t begin, end;
    begin = clock();
    for (int i = 0; i < n; i++) {
        EllipticCurve::dbl(R, G);
    }
    end = clock();
    printf("\tJacobi dbl\ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}

void benchmark() {
    mp_limb_t p[SIZE] = {
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0x1ff
    };

    mp_limb_t a[SIZE] = {
        0x2137201477420138L, 0x6313232130472104L, 0x3639164821638291L, 0x8371372917431694L,
        0x1739278732747374L, 0x1732913789724789L, 0x7296392641693692L, 0x32dead13L,
    };

    mp_limb_t b[SIZE] = {
        0x8032810382412080L, 0x2427104721732301L, 0x3281037402174017L, 0x7201839284217487L, 
        0x1737201372138214L, 0x3747217427103782L, 0x13729dee74L,
    };
    mp_limb_t c[SIZE] = {0};

    mp_limb_t tmp_z[SIZE * 2] = {0};
    mp_limb_t t[SIZE*2] = {0};
    mp_limb_t s[SIZE*2] = {0};


    const int n = 100000;
    time_t begin, end;
    begin = clock();
    for (int i = 0; i < n; i++) {
        mpn_mul_n(tmp_z, (const mp_limb_t *)a, (const mp_limb_t *)b, SIZE);
        mod((mp_limb_t*)c, (const mp_limb_t *)tmp_z, (const mp_limb_t *)p, t, s);
    }
    end = clock();
    printf("\tmulMod \ttime = %fusec\n", (end - begin) / double(CLOCKS_PER_SEC) / n * 1e6);
}


int main() {
#ifdef SECP521
    benchmark();
#endif
    benchmark_FP521mul();
    benchmark_FP521sqr();
    benchmark_secp521r_dbl();
    benchmark_secp521r();
}
