

class Fp {
public:
    static mpz_class modulus;
    mpz_class value; // 0 <= value < modulus

    Fp() { }
    Fp(mpz_class v) : value(v){ }

    static void setModulo(mpz_class v) {
        modulus = v;
    }

    Fp operator+(const Fp& other) {
        Fp r;
        r.value = (value + other.value) % modulus;
        if (r.value < 0) {
            r.value += modulus;
        }
        return r;
    }

    Fp operator-(const Fp& other) {
        Fp r;
        r.value = (value - other.value) % modulus;
        if (r.value < 0) {
            r.value += modulus;
        }
        return r;
    }

    Fp operator*(const Fp& other) {
        Fp r;
        r.value = (value * other.value) % modulus;
        if (r.value < 0) {
            r.value += modulus;
        }
        return r;
    }
};

mpz_class Fp::modulus;
