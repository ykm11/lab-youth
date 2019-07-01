#include<stdio.h>
#include<gmpxx.h>
#include<iostream>

#include "FP.h"

class Point;
void add(Point& R, const Point& P, const Point& Q);

class Point {
public:
    Fp x, y, z;

    Point() {}
    Point(const Fp& X, const Fp& Y, const Fp& Z) : x(X), y(Y), z(Z) {}


    Point operator+(const Point& other) const {
        Point r;
        add(r, *this, other);
        return r;
    }
};


class EllipticCurve {
public:
    Fp a, b;

    EllipticCurve() {}
    EllipticCurve(const mpz_class& A, const mpz_class& B) {
        a = Fp(A);
        b = Fp(B);
    }

    Point point(const mpz_class& x, const mpz_class& y) const {
        Fp r, l, u;
        mul(l, y, y); // y^2
        mul(r, x, x);
        mul(r, r, x); // x^3
        mul(u, a, x); // ax
        add(r, r, u); // x^3 + ax 
        add(r, r, b); // x^3 + ax + b

        Point P;
        if (l == r) {
            P = Point(x, y, one);
        } else { // 曲線に乗らない場合をどうするか
            P = Point(zero, one, zero);
        }
        return P;
    }
};

void add(Point& R, const Point& P, const Point& Q) {
    Fp u, v, v2, v3, w, s, t;
    Fp Rx, Ry, Rz;

    mul(s, Q.y, P.z); // Y2 * Z1
    mul(t, P.y, Q.z); // Y1 * Z2
    sub(u, s, t); // Y2 * Z1 - Y1 * Z2

    mul(s, Q.x, P.z); // X2 * Z1
    mul(t, P.x, Q.z); // X1 * Z2
    sub(v, s, t); // X2 * Z1 - X1 * Z2

    mul(v2, v, v); // v^2
    mul(v3, v2, v); // v^3

    mul(t, two, v2); // 2 * v^2
    mul(t, t, P.x); // 2 * v^2 * X1
    mul(t, t, Q.z); //  2 * v^2 * X1*Z2

    mul(s, u, u); // u^2
    mul(s, s, P.z); // u^2 * Z1
    mul(s, s, Q.z); // u^2 * Z1 * Z2
    
    sub(w, s, t); //  (u^2 * Z1 * Z2) - (2 * v^2 * X1*Z2)
    sub(w, w, v3); // (u^2 * Z1 * Z2) - v^3 - (2 * v^2 * X1*Z2)

    mul(Rx, v, w); // Rx = v * w

    mul(t, v3, P.y); // v^3 * Y1
    mul(t, t, Q.z); // v^3 * Y1 * Z2
    mul(s, v2, P.x); // v^2 * X1
    mul(s, s, Q.z); // v^2 * X1 * Z2
    sub(s, s, w); // v^2 * X1 * Z2 - w
    mul(s, s, u); // u(v^2 * X1 * Z2 - w)
    sub(Ry, s, t); // Ry = u(v^2 * X1 * Z2 - w) -  v^3 * Y1 * Z2

    mul(Rz, v3, P.z);
    mul(Rz, Rz, Q.z); // Rz = v^3 * Z1 * Z2

    if (Rz == zero) {
        R.x = zero;    
        R.y = one;    
        R.z = zero;    
    } else {
        R.x = Rx;
        R.y = Ry;
        R.z = Rz;
    }
}
