#pragma once
#include<stdio.h>
#include<gmpxx.h>
#include<iostream>

#include "FP.h"

class Point;
class EllipticCurve;

bool isEqual(const Point &P, const Point &Q);
void add(Point &R, const Point &P, const Point &Q);
void sub(Point &R, const Point &P, const Point &Q);
void mul(Point &R, const Point &P, const mpz_class &x); // scalar multiplying
void dbl(Point &R, const Point &P);

void dump(const Point &P);

void r_mul(Point &R, const Point &G, const mpz_class &x);
void montgomery_mul(Point &R0, const Point &G, const mpz_class &n);
void window_mul(Point &R, const Point &G, const mpz_class &n);

void multipleMul(Point &R, const Point &P, const mpz_class &u, const Point &Q, const mpz_class &v);

void naf_mul(Point &R, const Point &G, const mpz_class &x);

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
        mul(r, *this, x);
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


class EllipticCurve {
public:
    static Fp a, b;

    EllipticCurve() {}
    EllipticCurve(const mpz_class &A, const mpz_class &B) {
        a = Fp(A);
        b = Fp(B);
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
};


class GLV {

public:
    static Fp rw; 
    static mpz_class lmd, order; 
    static Point base; 

    static void initForsecp256k1() {
        Fp::setModulo(mpz_class("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16));
        Fp::squareRoot(rw, Fp(-3));
        mpz_add_ui(rw.value.get_mpz_t(), rw.value.get_mpz_t(), 1);
        Fp::neg(rw, rw); // - (sqrt(-3) + 1)
        mpz_tdiv_q_2exp(rw.value.get_mpz_t(), rw.value.get_mpz_t(), 1); 
        // - (sqrt(-3) + 1) / 2

        Fp gx, gy;
        gx = Fp("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
        gy = Fp("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
        base = Point(gx.value, gy.value, 1);
        
        order = mpz_class("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16);
        lmd = mpz_class("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16);

    }

    static void decomposing_k(mpz_class &k0, mpz_class &k1, const mpz_class &k);
    static void mulBase(Point &R, const mpz_class &k);
    static void lambdaMul(Point &R, const Point &P);
    static void scalarMul(Point &R, const Point &P, const mpz_class &k);
};


static inline void getNafArray(int8_t naf[], const mpz_class &x) {
    int j = 0;
    mpz_class z, n; 
    n = x;
    while (n > 0) {
        if((n & 1) == 1) {
            z = 2 - (n & 3);
            n = n - z;
        } else {
            z = 0;
        }
        n = n >> 1;
        naf[j] = mpz_get_si(z.get_mpz_t());
        j++;
    }
}
