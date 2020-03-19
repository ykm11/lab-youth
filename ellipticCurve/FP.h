#pragma once
#include<gmpxx.h>

#include <iostream>

class Fp;
void add(Fp& z, const Fp& x, const Fp& y);
void sub(Fp& z, const Fp& x, const Fp& y);
void mul(Fp& z, const Fp& x, const Fp& y);
void invmod(Fp& r, const Fp& x);
bool isEq(const Fp& x, const Fp& y);
void sqr(Fp& r, const Fp& x);

static inline bool zeroCmp(const Fp &x);

#ifdef USE_MPN
#define SIZE 4

void add(Fp& z, const Fp& x, uint64_t scalar);
#endif

class Fp {
public:
#ifndef USE_MPN
    static mpz_class modulus;
    mpz_class value; // 0 <= value < modulus

    Fp() { }
    Fp(const mpz_class& v) : value(v % modulus) {
        if(value < 0) {
            value += modulus;
        }
    }

    Fp(const std::string& str, int base) : value(mpz_class(str, base) % modulus) {
        if(value < 0) {
            value += modulus;
        }
    }

    static void setModulo(const mpz_class& v);
#else
    static uint64_t modulus[SIZE];
    uint64_t value[SIZE] = {0};

    Fp() { }
    Fp(uint64_t v[SIZE]){
        for (size_t i = 0; i < SIZE; ++i) {
            value[i] = v[i];
        }
        if (mpn_cmp((const mp_limb_t *)value, (const mp_limb_t *)modulus, SIZE) >= 0) {
            mpn_sub_n((mp_limb_t *)value, (const mp_limb_t *)value, (const mp_limb_t *)Fp::modulus, SIZE);
        }
    }

    Fp(const mpz_class& v) {
        mpz_class x = v;
        for (size_t i = 0; i < SIZE; i++) {
            value[i] = mpz_get_ui(x.get_mpz_t());
            mpz_tdiv_q_2exp(x.get_mpz_t(), x.get_mpz_t(), 64);
        }
        if (mpn_cmp((const mp_limb_t *)value, (const mp_limb_t *)modulus, SIZE) >= 0) {
            mpn_sub_n((mp_limb_t *)value, (const mp_limb_t *)value, (const mp_limb_t *)modulus, SIZE);
        }
    }

    static void setModulo(const uint64_t p[SIZE]);
    static void setModulo(const mpz_class& v);

#endif

    Fp operator+(const Fp& other) const { 
        Fp z; 
        add(z, *this, other); 
        return z; 
    }

    Fp operator-(const Fp& other) const {
        Fp z; 
        sub(z, *this, other); 
        return z; 
    }

    Fp operator*(const Fp& other) const {
        Fp z; 
        mul(z, *this, other); 
        return z; 
    }

    bool operator==(const Fp& other) const {
        return isEq(*this, other);
    }
    bool operator!=(const Fp& other) const {
        return !isEq(*this, other);
    }

    static void neg(Fp& r, const Fp& x);
    static void mulInt(Fp& z, const Fp& x, int scalar);
    static bool squareRoot(Fp& r, const Fp& x);

};


static inline void move(Fp &z, const Fp &x) {
#ifndef USE_MPN
    z.value = x.value;
#else
    for (size_t i = 0; i < SIZE; i++) {
        z.value[i] = x.value[i];
    }
#endif
}


static inline bool zeroCmp(const Fp &x) {
#ifndef USE_MPN
    return (x.value == 0);
#else
    return (mpn_zero_p((const mp_limb_t *)x.value, SIZE) == 1);
#endif
}


#ifdef USE_MPN
static inline int cmp(const uint64_t x[SIZE], const uint64_t y[SIZE]) {
    return mpn_cmp((const mp_limb_t *)x, (const mp_limb_t *)y, SIZE);
}


static inline void set_mpz_t(mpz_t& z, const uint64_t* p, int n) {
    int s = n;
    while (s > 0) {
        if (p[s - 1]) break;
        s--;
    }
    z->_mp_alloc = n;
    z->_mp_size = s;
    z->_mp_d = (mp_limb_t*)const_cast<uint64_t*>(p);
}

static inline void dump(const mp_limb_t x[SIZE]) {
    mpz_t mx;
    set_mpz_t(mx, (const uint64_t*)x, SIZE+1);
    std::cout << mx << std::endl;
}

#endif

static inline void dump(const Fp &x) {
#ifdef USE_MPN
    mpz_t mx;
    set_mpz_t(mx, x.value, SIZE);
    std::cout << mx << std::endl;
#else
    std::cout << x.value << std::endl;
#endif
}



static inline void mulMod(mpz_class& z, const mpz_class& x, const mpz_class& y, const mpz_class& m) {
    mpz_mul(z.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
    mpz_mod(z.get_mpz_t(), z.get_mpz_t(), m.get_mpz_t());
}

static inline void sqrMod(mpz_class& z, const mpz_class& x, const mpz_class& m) {
    mpz_powm_ui(z.get_mpz_t(), x.get_mpz_t(), 2, m.get_mpz_t());
}

static inline void powMod(mpz_class& z, const mpz_class& x, const mpz_class& y, const mpz_class& m) {
    mpz_powm(z.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), m.get_mpz_t());
}


