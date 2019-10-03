#pragma once
#include<gmpxx.h>

class Fp;
void add(Fp& z, const Fp& x, const Fp& y);
void sub(Fp& z, const Fp& x, const Fp& y);
void mul(Fp& z, const Fp& x, const Fp& y);
void mul(Fp& z, const Fp& x, int scalar);
void invmod(Fp& r, const Fp& x);
bool isEq(const Fp& x, const Fp& y);
void sqr(Fp& r, const Fp& x);


class Fp {
public:
    static mpz_class modulus;
    mpz_class value; // 0 <= value < modulus

    Fp() { }
#if 1
    Fp(const mpz_class& v) : value(v % modulus) {
        if(value < 0) {
            value += modulus;
        }
    }
#else
    Fp(const mpz_class& v) : value(v) {
        if(v >= modulus) {
            value %= modulus;
        }
        if(value < 0) {
            value %= modulus;
            value += modulus;
        }
    }
#endif
    //Fp(mpz_class v) : value(std::move(v)) { }

    Fp(const std::string& str, int base) : value(mpz_class(str, base) % modulus) { // intにconst はつけなくてよい
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

    static void mulInt(Fp& z, const Fp& x, int scalar);

};


