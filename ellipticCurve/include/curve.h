#pragma once
#include<stdio.h>
#include<gmpxx.h>
#include<iostream>

#include "FP.h"

class Point;
class jPoint;
class EllipticCurve;

bool isEqual(const Point&, const Point&);
bool isEqual(const jPoint&, const jPoint&);
void add(Point&, const Point&, const Point&);
void add(jPoint&, const jPoint&, const jPoint&);
template<class TPoint> void sub(TPoint&, const TPoint&, const TPoint&);

template<class TPoint> void dump(const TPoint P);

template<class TPoint> void l_mul(TPoint&, const TPoint&, const mpz_class&); 
template<class TPoint> void r_mul(TPoint&, const TPoint&, const mpz_class&);
template<class TPoint> void montgomery_mul(TPoint&, const TPoint&, const mpz_class&);
template<class TPoint> void window_mul(TPoint&, const TPoint&, const mpz_class&);
template<class TPoint> void naf_mul(TPoint&, const TPoint&, const mpz_class&);

void multipleMul(Point&, const Point&, const mpz_class&, const Point&, const mpz_class&);

inline void getNafArray(int8_t *naf, const mpz_class&);
template <class TPoint> inline void setPoint(TPoint&, const Fp&, const Fp&, const Fp&);

class Point {
public:
    Fp x, y, z;

    Point() {}
    Point(const Fp &X, const Fp &Y, const Fp &Z) : x(X), y(Y), z(Z) {}
    Point(const mpz_class &X, const mpz_class &Y, const mpz_class &Z) : x(X), y(Y), z(Z) {}
    //~Point() = default;


    Point operator+(const Point &other) const {
        Point r;
        add(r, *this, other);
        return r;
    }

    Point operator-(const Point &other) const {
        Point r;
        sub(r, *this, other);
        return r;
    }

    Point operator*(const mpz_class &x) const {
        Point r;
        naf_mul(r, *this, x);
        return r;
    }

    bool operator==(const Point &other) const {
        return isEqual(*this, other);
    }

    bool operator!=(const Point &other) const {
        return !isEqual(*this, other);
    }

    void xy(Fp &s, Fp &t) const { // 射影座標からアフィン座標へ
        Fp inv_z;
        invmod(inv_z, z);
        
        mul(s, x, inv_z); // s <- X/Z
        mul(t, y, inv_z); // t <- Y/Z
    }

    static void neg(Point &R, const Point &P) {
        R.x = P.x;
        R.z = P.z;
        Fp::neg(R.y, P.y); 
    }
};


class jPoint { // Jacobian coordinate
public:
    Fp x, y, z;

    jPoint() {}
    jPoint(const Fp &X, const Fp &Y, const Fp &Z) : x(X), y(Y), z(Z) {}
    jPoint(const mpz_class &X, const mpz_class &Y, const mpz_class &Z) : x(X), y(Y), z(Z) {}

    void xy(Fp &s, Fp &t) const { // ヤコビ座標からアフィン座標へ
        Fp inv_z;
        invmod(inv_z, z);

        sqr(t, inv_z); // Z^{2}
        mul(s, x, t); // X / Z^{2}
        mul(t, t, inv_z); // Z^{3}
        mul(t, y, t); // Y / Z^{3}
    }

    jPoint operator+(const jPoint &other) const {
        jPoint r;
        add(r, *this, other);
        return r;
    }

    jPoint operator-(const jPoint &other) const {
        jPoint r;
        sub(r, *this, other);
        return r;
    }

    jPoint operator*(const mpz_class &x) const {
        jPoint r;
        naf_mul(r, *this, x);
        return r;
    }

    bool operator==(const jPoint &other) const {
        return isEqual(*this, other);
    }

    bool operator!=(const jPoint &other) const {
        return !isEqual(*this, other);
    }

    static void neg(jPoint &R, const jPoint &P) {
        R.x = P.x;
        R.z = P.z;
        Fp::neg(R.y, P.y); 
    }
};

class EllipticCurve {
public:
    static Fp a, b;

    EllipticCurve() {}
    EllipticCurve(const mpz_class &A, const mpz_class &B) {
        a = Fp(A);
        b = Fp(B);
    }
    EllipticCurve(const Fp &A, const Fp &B) {
        copy(a, A);
        copy(b, B);
    }
    //~EllipticCurve() = default;

    Point point(const mpz_class &x, const mpz_class &y) const {
        Point P = point(x, y, 1);
        return P;
    }

    Point point(const mpz_class &x, const mpz_class &y, const mpz_class &z) const {
        Point P;
        Fp r, l, u, v, z2;

        sqr(z2, z);

        sqr(l, y);
        mul(l, l, z); // Z * Y^2
        sqr(r, x);
        mul(r, r, x); // X^3
        mul(u, a, x); // aX
        mul(u, u, z2); // aX * Z^2
        add(r, r, u); // X^3 + aX * Z^2
        mul(v, b, z);
        mul(v, v, z2); // b * Z^3
        add(r, r, v); // X^3 + aX * Z^2 + b * Z^3

        if (l == r) {
            P = Point(x, y, z);
        } else { 
            throw "Exception: Point does not exist on the curve";
        }
        return P;
    }
 
    Point operator()(const mpz_class &x, const mpz_class &y) const {
        Point P = point(x, y, 1);
        return P;
    }

    Point operator()(const mpz_class &x, const mpz_class &y, const mpz_class &z) const {
        Point P = point(x, y, z);
        return P;
    }

    static void dbl(Point &R, const Point &P);
    static void dbl(jPoint &R, const jPoint &P);
};

template<class TPoint>
void sub(TPoint &R, const TPoint &P, const TPoint &Q) {
    TPoint::neg(R, Q); // R <- [-1]Q
    add(R, P, R); // R <- P + [-1]Q
}

template<class TPoint>
void l_mul(TPoint &R, const TPoint &P, const mpz_class &x) { // 左向きバイナリ法
    TPoint tmp_P = P;
    setInfPoint(R);

    size_t k_bits = mpz_sizeinbase(x.get_mpz_t(), 2)-1;
    for (size_t i = 0; i < k_bits; i++) {
        if ((mpz_tstbit(x.get_mpz_t(), i)) == 1) {
            add(R, R, tmp_P);
        }
        EllipticCurve::dbl(tmp_P, tmp_P);
    }
    add(R, R, tmp_P);
}

template<class TPoint>
void r_mul(TPoint &R, const TPoint& G, const mpz_class &x) { // 右向きバイナリ法
    size_t k_bits = mpz_sizeinbase(x.get_mpz_t(), 2);
    setPoint(R, G.x, G.y, G.z);

    for (int i = k_bits-2; i >= 0; i--) {
        EllipticCurve::dbl(R, R);

        if (mpz_tstbit(x.get_mpz_t(), i) == 1) {
            add(R, R, G);
        }
    }
}

template<class TPoint>
void montgomery_mul(TPoint &R0, const TPoint &G, const mpz_class &n) { // Montgomery Ladder
    size_t k_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    TPoint R1 = G;

    setInfPoint(R0);
    for (int i = k_bits-1; i >= 0; i--) {
        if (mpz_tstbit(n.get_mpz_t(), i) == 0) {
            add(R1, R1, R0);
            EllipticCurve::dbl(R0, R0);
        } else {
            add(R0, R0, R1);
            EllipticCurve::dbl(R1, R1);
        }
    }
}

template<class TPoint>
void window_mul(TPoint &R, const TPoint &G, const mpz_class &n) { // windows method
    size_t k_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    TPoint P[4];

    setInfPoint(P[0]);
    P[1] = G;
    EllipticCurve::dbl(P[2], G);
    add(P[3], P[2], G);

    R = P[0];
    for (int i = k_bits-1; i > 0; i=i-2) {
        EllipticCurve::dbl(R, R);
        EllipticCurve::dbl(R, R); // R <- 4R

        add(R, R, P[2*mpz_tstbit(n.get_mpz_t(), i) + mpz_tstbit(n.get_mpz_t(), i-1)]);
    } 
    if ((k_bits & 1) == 1) { // nのビット数が奇数のときだけ
        EllipticCurve::dbl(R, R);
        add(R, R, P[mpz_tstbit(n.get_mpz_t(), 0)]);
    }
}

template<class TPoint>
void naf_mul(TPoint &R, const TPoint &P, const mpz_class &x) {
    size_t naf_size = mpz_sizeinbase(x.get_mpz_t(), 2) + 1;
    int8_t naf[naf_size];
    memset(naf, 0, naf_size);

    size_t w_size = 5;
    size_t tblSize = 1 << w_size;
    TPoint tbl[tblSize];
    setInfPoint(tbl[0]);
    tbl[1] = P;

    for (size_t k = 2; k < 21; k=k+2) {
        EllipticCurve::dbl(tbl[k], tbl[k/2]);
        add(tbl[k+1], tbl[k], P);
    }
    
    getNafArray(naf, x);
    while (naf_size >= 1 && naf[naf_size-1] == 0) {
        naf_size--;
    }

    TPoint Q;
    R = P;
    int8_t t;
    int i;
    for (i = naf_size-2; i >= int(w_size-1); i=i-w_size) {
        for (size_t j = 0; j < w_size; j++) {
            EllipticCurve::dbl(R, R);
        }

        t = 16*naf[i] + 8*naf[i-1] + 4*naf[i-2] + 2*naf[i-3] + naf[i-4]; 
        if (t < 0) {
            TPoint::neg(Q, tbl[-t]);
        } else {
            Q = tbl[t];
        }
        add(R, R, Q);
    }

    if (i < 0) {
        return;
    }

    t = 0;
    while (i > 0) {
        EllipticCurve::dbl(R, R);
        t = t + (1 << i) * naf[i];
        i--;
    }
    EllipticCurve::dbl(R, R);
    t = t + naf[0];

    if (t < 0) {
        TPoint::neg(Q, tbl[-t]);
    } else {
        Q = tbl[t];
    }
    add(R, R, Q);
}


class GLV {

public:
    static Fp rw; 
    static mpz_class lmd, order; 
    static Point base; 

    static void initForsecp256k1() {
        Fp::setModulo(mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16));
        Fp::squareRoot(rw, Fp(-3));
#ifndef YKM_ECC_USE_MPN
        mpz_add_ui(rw.value.get_mpz_t(), rw.value.get_mpz_t(), 1);
        Fp::neg(rw, rw); // - (sqrt(-3) + 1)
        mpz_tdiv_q_2exp(rw.value.get_mpz_t(), rw.value.get_mpz_t(), 1); 
        // - (sqrt(-3) + 1) / 2

        Fp gx(mpz_class("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16));
        Fp gy(mpz_class("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16));
        base = Point(gx.value, gy.value, 1);
#else
        add(rw, rw, 1);
        Fp::neg(rw, rw); // - (sqrt(-3) + 1)
        mpn_rshift((mp_limb_t *)rw.value, (const mp_limb_t *)rw.value, Fp::size_, 1);
        // - (sqrt(-3) + 1) / 2

        mpz_class gx("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
        mpz_class gy("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
        base = Point(gx, gy, 1);
#endif
        order = mpz_class("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16);
        lmd = mpz_class("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16);
    }

    static void decomposing_k(mpz_class &k0, mpz_class &k1, const mpz_class &k);
    static void mulBase(Point &R, const mpz_class &k);
    static void lambdaMul(Point &R, const Point &P);
    static void scalarMul(Point &R, const Point &P, const mpz_class &k);
};

template<class TPoint> void dump(const TPoint P) {
    if (zeroCmp(P.z)) {
        std::cout << "(0 : 1 : 0)" << std::endl;
    } else {
        Fp x, y;
        P.xy(x, y);
#ifndef YKM_ECC_USE_MPN
        std::cout << "(" << x.value << " : " << y.value << " : 1)" << std::endl;
#else
        mpz_t mx, my;
        set_mpz_t(mx, (uint64_t*)x.value, Fp::size_);
        set_mpz_t(my, (uint64_t*)y.value, Fp::size_);
        std::cout << "(" << mx << " : " << my << " : 1)" << std::endl;
#endif
    }
}

inline void setInfPoint(Point &R) {
#ifndef YKM_ECC_USE_MPN
    R.x.value = 0;
    R.y.value = 1;
    R.z.value = 0;
#else
    mpn_zero((mp_limb_t *)R.x.value, Fp::size_);
    mpn_zero((mp_limb_t *)R.y.value, Fp::size_);
    mpn_zero((mp_limb_t *)R.z.value, Fp::size_);
    //R.x.value[0] = 0;
    R.y.value[0] = 1;
    //R.z.value[0] = 0;
#endif
}

inline void setInfPoint(jPoint &R) {
#ifndef YKM_ECC_USE_MPN
    R.x.value = 1;
    R.y.value = 1;
    R.z.value = 0;
#else
    mpn_zero(R.x.value, Fp::size_);
    mpn_zero(R.y.value, Fp::size_);
    mpn_zero(R.z.value, Fp::size_);
    R.x.value[0] = 1;
    R.y.value[0] = 1;
    //R.z.value[0] = 0;
#endif
}

template <class TPoint> inline void setPoint(TPoint &R, const Fp &Rx, const Fp &Ry, const Fp &Rz) {
#ifndef YKM_ECC_USE_MPN
    R.x.value = Rx.value;
    R.y.value = Ry.value;
    R.z.value = Rz.value;
#else
    copy_n(R.x.value, (mp_limb_t *)Rx.value, Fp::size_);
    copy_n(R.y.value, (mp_limb_t *)Ry.value, Fp::size_);
    copy_n(R.z.value, (mp_limb_t *)Rz.value, Fp::size_);
#endif
}

inline void getNafArray(int8_t *naf, const mpz_class &x) {
    int j = 0;
    mpz_class z, n; 
    n = x;
    while (n > 0) {
        if((n & 1) == 1) {
            z = 2 - (n & 3);
            n -= z;
        } else {
            z = 0;
        }
        naf[j] = mpz_get_si(z.get_mpz_t());
        n >>= 1;
        j++;
    }
}

