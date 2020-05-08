#pragma once
#include<gmpxx.h>

#include <iostream>

class Fp;
void add(Fp&, const Fp&, const Fp&);
void sub(Fp&, const Fp&, const Fp&);
void mul(Fp&, const Fp&, const Fp&);
void invmod(Fp&, const Fp&);
bool isEq(const Fp&, const Fp&);
void sqr(Fp&, const Fp&);

inline void copy(Fp &z, const Fp &x);
inline bool zeroCmp(const Fp&);

#ifdef YKM_ECC_USE_MPN
#define YKM_ECC_MAX_SIZE ((521+63)/64)

inline mp_limb_t sub_n(mp_limb_t*, mp_limb_t*, mp_limb_t*, size_t);
inline mp_limb_t add_n(mp_limb_t*, mp_limb_t*, mp_limb_t*, size_t);
inline void mul_n(mp_limb_t*, mp_limb_t*, mp_limb_t*, size_t);
inline void getArray(mp_limb_t*, size_t, const mpz_class&, int);
inline void mulMod(mp_limb_t*, mp_limb_t*, mp_limb_t*, mp_limb_t*, 
        mp_limb_t*,  mp_limb_t*, size_t);
inline void powMod(mp_limb_t*, mp_limb_t*, mp_limb_t*, 
        mp_limb_t*, mp_limb_t*, size_t);
inline void sqrMod(mp_limb_t*, mp_limb_t*, mp_limb_t*, 
        mp_limb_t*,  mp_limb_t*, size_t);

inline int cmp_n(mp_limb_t *x, mp_limb_t *y, size_t n);
inline void copy_n(mp_limb_t *r, mp_limb_t *v, size_t n);

#else
inline void mulMod(mpz_class&, const mpz_class&, const mpz_class&, const mpz_class&);
inline void sqrMod(mpz_class&, const mpz_class&, const mpz_class&);
inline void powMod(mpz_class&, const mpz_class&, const mpz_class&, const mpz_class&);
#endif

 
#ifdef YKM_ECC_SECP521
void add(Fp&, const Fp&, uint64_t);
inline void secp521Mod(mp_limb_t*, const mp_limb_t*, const mp_limb_t*, mp_limb_t*, mp_limb_t*);
 
#elif defined(YKM_ECC_USE_MPN)
void add(Fp&, const Fp&, uint64_t);

#endif

#ifndef YKM_ECC_USE_MPN
class Fp {
public:
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
#else
class Fp {
public:
    static mp_limb_t modulus[YKM_ECC_MAX_SIZE];
    static size_t size;
    mp_limb_t value[YKM_ECC_MAX_SIZE];

    Fp() { }
    Fp(mp_limb_t v[YKM_ECC_MAX_SIZE]){
        mpn_copyi(value, (const mp_limb_t *)v, size);
        if (mpn_cmp((const mp_limb_t *)value, (const mp_limb_t *)modulus, size) >= 0) {
            sub_n(value, value, modulus, size);
        }
    }

    Fp(const mpz_class& v) {
        getArray(value, size, v, v.get_mpz_t()->_mp_size);

        if (v < 0) {
            sub_n(value, modulus, value, size);
            return;
        }
        if (mpn_cmp((const mp_limb_t *)value, (const mp_limb_t *)modulus, size) >= 0) {
            sub_n(value, value, modulus, size);
            return;
        }
    }

    static void setModulo(const mp_limb_t p[YKM_ECC_MAX_SIZE]);
    static void setModulo(const mpz_class& v);

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
#endif


inline void copy(Fp &z, const Fp &x) {
#ifndef YKM_ECC_USE_MPN
    z.value = x.value;
#else
    copy_n(z.value, (mp_limb_t *)x.value, Fp::size);
#endif
}


inline bool zeroCmp(const Fp &x) {
#ifndef YKM_ECC_USE_MPN
    return (x.value == 0);
#else
    return (mpn_zero_p(x.value, Fp::size) == 1);
#endif
}


#ifdef YKM_ECC_USE_MPN

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

#endif

static inline void dump(const Fp &x) {
#ifdef YKM_ECC_USE_MPN
    mpz_t mx;
    set_mpz_t(mx, (const uint64_t*)x.value, Fp::size);
    std::cout << mx << std::endl;
#else
    std::cout << x.value << std::endl;
#endif
}

#ifdef YKM_ECC_USE_MPN
inline void powMod(mp_limb_t* r, mp_limb_t* x, mp_limb_t* e, mp_limb_t* modulus,
        mp_limb_t* tp, size_t size) {
    mpn_sec_powm(r, (const mp_limb_t *)x, size, (const mp_limb_t*)e, 
            size*GMP_NUMB_BITS, (const mp_limb_t *)modulus, size, tp);
}

inline void sqrMod(mp_limb_t* r, mp_limb_t* x, mp_limb_t* modulus, 
        mp_limb_t* tmp,  mp_limb_t* q, size_t size) {
    mpn_sqr(tmp, (const mp_limb_t *)x, size);
    mpn_tdiv_qr(q, r, 0, (const mp_limb_t*)tmp, size*2, (const mp_limb_t *)modulus, size);
}

inline void mulMod(mp_limb_t* z, mp_limb_t* x, mp_limb_t* y, mp_limb_t* modulus, 
        mp_limb_t* tmp,  mp_limb_t* q, size_t size) {
    mul_n(tmp, x, y, size);
    mpn_tdiv_qr(q, z, 0, (const mp_limb_t*)tmp, size*2, (const mp_limb_t *)modulus, size);
} 

#else
inline void mulMod(mpz_class& z, const mpz_class& x, const mpz_class& y, const mpz_class& m) {
    mpz_mul(z.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
    mpz_mod(z.get_mpz_t(), z.get_mpz_t(), m.get_mpz_t());
}

inline void sqrMod(mpz_class& z, const mpz_class& x, const mpz_class& m) {
    mpz_powm_ui(z.get_mpz_t(), x.get_mpz_t(), 2, m.get_mpz_t());
}

inline void powMod(mpz_class& z, const mpz_class& x, const mpz_class& y, const mpz_class& m) {
    mpz_powm(z.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), m.get_mpz_t());
}

#endif

inline mp_limb_t add_n(mp_limb_t* z, mp_limb_t* x, mp_limb_t* y, size_t n) {
    return mpn_add_n(z, (const mp_limb_t *)x, (const mp_limb_t*)y, n); 
}

inline mp_limb_t sub_n(mp_limb_t* z, mp_limb_t* x, mp_limb_t* y, size_t n) {
    return mpn_sub_n(z, (const mp_limb_t *)x, (const mp_limb_t*)y, n); 
}

inline void mul_n(mp_limb_t* XY, mp_limb_t* x, mp_limb_t* y, size_t n) {
    mpn_mul_n(XY, (const mp_limb_t*)x, (const mp_limb_t*)y, n);
}

inline void copy_n(mp_limb_t *r, mp_limb_t *v, size_t n) {
    mpn_copyi(r, (const mp_limb_t*)v, n);
}

inline int cmp_n(mp_limb_t *x, mp_limb_t *y, size_t n) {
    return mpn_cmp((const mp_limb_t*)x, (const mp_limb_t*)y, n);
}


inline void getArray(mp_limb_t *buf, size_t maxSize, const mpz_class &x, int xn) {
    const size_t bufByteSize = sizeof(mp_limb_t) * maxSize;
    size_t xByteSize;
    if (xn < 0) {
        xByteSize = sizeof(mp_limb_t) * (-xn);
    } else {
        xByteSize = sizeof(mp_limb_t) * xn;
    }
    memcpy(buf, x.get_mpz_t()->_mp_d, xByteSize);
    memset((char*)buf + xByteSize, 0, bufByteSize - xByteSize);
}

#ifdef YKM_ECC_SECP521
inline void secp521Mod(mp_limb_t *z, const mp_limb_t *XY, const mp_limb_t *p, mp_limb_t *t, mp_limb_t *s) {
    // (T + (T mod R)*N) / R
    mpn_zero(t, 9*2);
    mpn_zero(s, 9*2);
    mpn_and_n(t, XY, p, 9); // T mod R
    for (size_t i = 0; i < 9; i++) {
        s[i+8] = t[i];
    }
    mpn_lshift(s, (const mp_limb_t*)s, 9*2, 9);
    sub_n(s, s, t, 9*2);
    add_n(t, s, (mp_limb_t *)XY, 9*2); // (T + (T mod R)*N)

    mpn_rshift(t, (const mp_limb_t*)t, 9*2, 9);
    for (size_t i = 0; i < 9; i++) { // (T + (T mod R)*N) / R
        t[i] = t[i+8];
    }
    if (mpn_cmp((const mp_limb_t*)t, p, 9) >= 0) {
        sub_n(t, t, (mp_limb_t*)p, 9);
    }
    copy_n(z, t, 9);
}
#endif
