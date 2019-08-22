#pragma once
#include<gmpxx.h>

class Fp;
void add(Fp& z, const Fp& x, const Fp& y);
void sub(Fp& z, const Fp& x, const Fp& y);
void mul(Fp& z, const Fp& x, const Fp& y);
void invmod(Fp& r, const Fp& x);
bool isEq(const Fp& x, const Fp& y);


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
    //Fp(mpz_class v) : value(std::move(v)) { }

    //Fp(const std::string& str, const int base) {
    Fp(const std::string& str, int base) { // const はつけなくてよい
        Fp(mpz_class(str, base));
    }

    static void setModulo(const mpz_class& v) {
        modulus = v;
    }

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

    static void mulInt(Fp& z, const Fp& x, int scalar) {
        z.value = x.value * scalar;
        z.value %= modulus;
    }

};
mpz_class Fp::modulus;

void add(Fp& z, const Fp& x, const Fp& y) {
    z.value = x.value + y.value;
    if(z.value >= Fp::modulus) {
        z.value -= Fp::modulus;
    }
}

void sub(Fp& z, const Fp& x, const Fp& y) {
    z.value = x.value - y.value;
    if (z.value < 0) {
        z.value += Fp::modulus;
    }
}

void mul(Fp& z, const Fp& x, const Fp& y) {
    z.value = (x.value * y.value) % Fp::modulus;
}

void invmod(Fp& r, const Fp& x) { 
    mpz_invert(r.value.get_mpz_t(), x.value.get_mpz_t(), Fp::modulus.get_mpz_t());
}

bool isEq(const Fp& x, const Fp& y) {
    return x.value == y.value;
}

