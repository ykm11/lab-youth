#include<stdio.h>
#include<gmpxx.h>
#include<iostream>

#include "FP.h"

class Point;
bool isEqual(const Point& P, const Point& Q);
void add(Point& R, const Point& P, const Point& Q);
void sub(Point& R, const Point& P, const Point& Q);
void mul(Point& R, const Point& P, const Fp& x); // scalar multiplying

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

    Point operator-(const Point& other) const {
        Point r;
        sub(r, *this, other);
        return r;
    }

    bool operator==(const Point& other) const {
        return isEqual(*this, other);
    }

    Point operator*(const Fp& x) const { // 右からかける P*x
        Point r;
        mul(r, *this, x);
        return r;
    }

    Point operator*(const mpz_class& x) const {
        Point r;
        Fp d = Fp(x);
        mul(r, *this, d);
        return r;
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
Fp EllipticCurve::a;
Fp EllipticCurve::b;


void add(Point& R, const Point& P, const Point& Q) {
    if (P.z == zero) {
        R.x = Q.x;
        R.y = Q.y;
        R.z = Q.z;
    } else if (Q.z == zero) {
        R.x = P.x;
        R.y = P.y;
        R.z = P.z;
    } else {
        Fp Rx, Ry, Rz;
        if (isEqual(P, Q)) { // if P == Q then Doubling
            Fp u, v, v2, w, s, t;

            mul(s, three, P.x); // 3*X
            mul(s, s, P.x); // 3*X^2
            mul(t, EllipticCurve::a, P.z); // a*Z
            mul(t, t, P.z); // a*Z^2
            add(u, s, t); // 3*X^2 + a*Z^2
        
            mul(v, P.y, P.z); // Y*Z

            mul(s, u, u); // u^2
            mul(t, eight, P.x); // 8*X
            mul(t, t, P.y); // 8*X*Y
            mul(t, t, v); // 8*X*Y*v
            sub(w, s, t); // u^2 - 8*X*Y*v

            mul(Rx, two, v); // 2*v
            mul(Rx, Rx, w); // Rx = 2*v*w

            mul(s, four, P.x);
            mul(s, s, P.y);
            mul(s, s, v);
            sub(s, s, w); // 4*X*Y*v - w
            mul(s, s, u); // u(4*X*Y*v - w)

            mul(v2, v, v); // v^2

            mul(t, P.y, P.y);
            mul(t, t, v2); // (Yv)^2
            mul(t, t, eight); // 8(Yv)^2

            sub(Ry, s, t); // Ry = u(4*X*Y*v - w) - 8(Yv)^2

            mul(Rz, eight, v); // 8*v
            mul(Rz, Rz, v2); // Rz = 8*v^3

        } else { // otherwise, Adding
            Fp u, v, v2, v3, w, s, t;

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
        }

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
}

void sub(Point& R, const Point& P, const Point& Q) {
    Point minus_Q = Q;
    sub(minus_Q.y, zero, minus_Q.y); // y <- p - y
    add(R, P, minus_Q);
}

void mul(Point& R, const Point& P, const Fp& x) {
    Point k = Point(zero, one, zero);
    Point tmp_P = P;

    Fp n = x;
    while (n.value > 0) {
        if (n.value % 2 == 1) {
            add(k, k, tmp_P);
        }
        
        add(tmp_P, tmp_P, P);
        n.value >>= 1;
    }
    R.x = k.x;
    R.y = k.y;
    R.z = k.z;
}

bool isEqual(const Point& P, const Point& Q) {
    Fp s, t, u, v;

    mul(s, Q.x, P.y);
    mul(t, P.x, Q.y);
    
    mul(u, P.x, Q.z);
    mul(v, Q.x, P.z);
 
    return (s == t) && (u == v);
}

