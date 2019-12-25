#pragma once
#include<stdio.h>
#include<gmpxx.h>
#include<iostream>

#include "FP.h"

class Point;
class EllipticCurve;

bool isEqual(const Point& P, const Point& Q);
void add(Point& R, const Point& P, const Point& Q);
void sub(Point& R, const Point& P, const Point& Q);
void mul(Point& R, const Point& P, const mpz_class& x); // scalar multiplying
void dbl(Point& R, const Point& P);

void print(const Point& P);

void r_mul(Point &R, const Point& G, const mpz_class x);
void montgomery_mul(Point &R0, const Point& G, const mpz_class n);
void window_mul(Point &R, const Point& G, const mpz_class n);

void multipleMul(Point &R, const Point &P, const mpz_class &u, const Point &Q, const mpz_class &v);


class Point {
public:
    Fp x, y, z;

    Point() {}
    Point(const Fp& X, const Fp& Y, const Fp& Z) : x(X), y(Y), z(Z) {}
    Point(const mpz_class& X, const mpz_class& Y, const mpz_class& Z) : x(X), y(Y), z(Z) {}
    //~Point() = default;


    Point operator+(const Point& other) const {
        Point r;
        add(r, *this, other);
        return r;
    }

    Point operator-(const Point& other) const {
        Point r;
        sub(r, *this, other);
        return r;
    }

    Point operator*(const mpz_class& x) const {
        Point r;
        mul(r, *this, x);
        return r;
    }

    bool operator==(const Point& other) const {
        return isEqual(*this, other);
    }

    bool operator!=(const Point& other) const {
        return !isEqual(*this, other);
    }

    void xy(Fp& s, Fp& t) const { // 射影座標からアフィン座標へ
        Fp inv_z;
        invmod(inv_z, z);
        
        mul(s, x, inv_z); // s <- X/Z
        mul(t, y, inv_z); // t <- Y/Z
    }
};


class EllipticCurve {
public:
    static Fp a, b;

    EllipticCurve() {}
    EllipticCurve(const mpz_class& A, const mpz_class& B) {
        a = Fp(A);
        b = Fp(B);
    }
    //~EllipticCurve() = default;

    Point point(const mpz_class& x, const mpz_class& y) const {
        Point P = point(x, y, 1);
        return P;
    }

    Point point(const mpz_class& x, const mpz_class& y, const mpz_class& z) const {
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
 
    Point operator()(const mpz_class& x, const mpz_class& y) const {
        Point P = point(x, y, 1);
        return P;
    }

    Point operator()(const mpz_class& x, const mpz_class& y, const mpz_class& z) const {
        Point P = point(x, y, z);
        return P;
    }

    static void dbl(Point& R, const Point& P);
};


class GLV {

public:
    static Fp rw; 
    static mpz_class lmd; 
    static Point base, base_; 

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
        mul(base_.x, rw, base.x);
        base_.y = base.y;
        base_.z = base.z;
        //print(base);
        //print(base_);
        
        lmd = mpz_class("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16);
        /*
            lmd = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72 
            [lmd]G = (rw * gx, gy) 
        */
    }

    static void mulBase(Point &R0, const mpz_class &k) { 
        // TODO p5 Algorithm 1
        mpz_class k0, k1; 
        /* k = k0 + k1 * lmd
           s.t 0 < k0, k1 < sqrt(|E|)
        */
        mpz_tdiv_qr(k1.get_mpz_t(), k0.get_mpz_t(), 
                k.get_mpz_t(), lmd.get_mpz_t()); 
        Point R1;
        mul(R1, base_, k1);
        mul(R0, base, k0);
        add(R0, R0, R1);
    }

};
