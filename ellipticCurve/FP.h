
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
        bool isEq;
        //isEq = (value == other.value);
        Fp z;
        sub(z, *this, other);
        isEq = (z.value == 0);
        return isEq;
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

Fp eight = Fp(8);
Fp four = Fp(4);
Fp three = Fp(3);
Fp two = Fp(2);
Fp one = Fp(1);
Fp zero = Fp(0);
