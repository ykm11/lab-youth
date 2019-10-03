#pragma once
#include<stdio.h>
#include<gmpxx.h>
#include<iostream>

#include "FP.h"

class Point;
bool isEqual(const Point& P, const Point& Q);
void add(Point& R, const Point& P, const Point& Q);
void sub(Point& R, const Point& P, const Point& Q);
void mul(Point& R, const Point& P, const mpz_class& x); // scalar multiplying
void dbl(Point& R, const Point& P);

void print(const Point& P);
void r_mul(Point &R, const Point& G, const mpz_class x);
void montgomery_mul(Point &R0, const Point& G, const mpz_class n);
void window_mul(Point &R, const Point& G, const mpz_class n);


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
        Point P;
#if 0
        Fp r, l, u;
        sqr(l, y); // y^2
        sqr(r, x);
        mul(r, r, x); // x^3
        mul(u, a, x); // ax
        add(r, r, u); // x^3 + ax 
        add(r, r, b); // x^3 + ax + b

        if (l == r) {
            P = Point(x, y, 1);
            //P = Point(x, y, Fp(1));
        } else { // 曲線に乗らない場合はエラーを返すか例外を投げる
            throw "Exception: Point does not exist on the curve";
        }
#else
        P = point(x, y, 1);
#endif
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

