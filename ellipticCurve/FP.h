
class Fp;
void add(Fp& z, const Fp& x, const Fp& y);
void sub(Fp& z, const Fp& x, const Fp& y);
void mul(Fp& z, const Fp& x, const Fp& y);

class Fp {
public:
    static mpz_class modulus;
    mpz_class value; // 0 <= value < modulus

    Fp() { }
    Fp(const mpz_class& v) : value(v) { }
    //Fp(mpz_class v) : value(std::move(v)) { }
     

    static void setModulo(mpz_class v) {
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
};
mpz_class Fp::modulus;

void add(Fp& z, const Fp& x, const Fp& y) {
    z.value = (x.value + y.value) % z.modulus;
    if (z.value < 0) {
        z.value += z.modulus;
    }
}

void sub(Fp& z, const Fp& x, const Fp& y) {
    z.value = (x.value - y.value) % z.modulus;
    if (z.value < 0) {
        z.value += z.modulus;
    }
}

void mul(Fp& z, const Fp& x, const Fp& y) {
    z.value = (x.value * y.value) % z.modulus;
    if (z.value < 0) {
        z.value += z.modulus;
    }
}


Fp two = Fp(2);
