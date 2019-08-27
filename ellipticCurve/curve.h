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

    static Point point(const mpz_class& x, const mpz_class& y) {
        Point P;
        Fp r, l, u;

        mul(l, y, y); // y^2
        mul(r, x, x);
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
        return P;
    }

    static Point point(const mpz_class& x, const mpz_class& y, const mpz_class& z) {
        Point P;
        Fp r, l, u, v, z2;

        mul(z2, z, z);

        mul(l, y, y);
        mul(l, l, z); // Z * Y^2
        mul(r, x, x);
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
        Point P = EllipticCurve::point(x, y);
        return P;
    }

    Point operator()(const mpz_class& x, const mpz_class& y, const mpz_class& z) const {
        Point P = EllipticCurve::point(x, y, z);
        return P;
    }


    static void dbl(Point& R, const Point& P) {
        if (P.z.value == 0) {
            R.x.value = 0;
            R.y.value = 1;
            R.z.value = 0;
        } else {
            Fp Rx, Ry, Rz;
            Fp u, v, v2, w, s, t;

            Fp::mulInt(s, P.x, 3);
            //mul(s, three, P.x); // 3*X
            mul(s, s, P.x); // 3*X^2
            mul(t, a, P.z); // a*Z
            mul(t, t, P.z); // a*Z^2
            add(u, s, t); // 3*X^2 + a*Z^2
        
            mul(v, P.y, P.z); // Y*Z

            mul(s, u, u); // u^2
            Fp::mulInt(t, P.x, 8);
            //mul(t, eight, P.x); // 8*X
            mul(t, t, P.y); // 8*X*Y
            mul(t, t, v); // 8*X*Y*v
            sub(w, s, t); // u^2 - 8*X*Y*v

            Fp::mulInt(Rx, v, 2);
            //mul(Rx, two, v); // 2*v
            mul(Rx, Rx, w); // Rx = 2*v*w

            Fp::mulInt(s, P.x, 4);
            //mul(s, four, P.x);
            mul(s, s, P.y);
            mul(s, s, v);
            sub(s, s, w); // 4*X*Y*v - w
            mul(s, s, u); // u(4*X*Y*v - w)

            mul(v2, v, v); // v^2

            mul(t, P.y, P.y);
            mul(t, t, v2); // (Yv)^2
            Fp::mulInt(t, t, 8);
            //mul(t, t, eight); // 8(Yv)^2

            sub(Ry, s, t); // Ry = u(4*X*Y*v - w) - 8(Yv)^2

            Fp::mulInt(Rz, v, 8);
            //mul(Rz, eight, v); // 8*v
            mul(Rz, Rz, v2); // Rz = 8*v^3

            if (Rz.value == 0) {
                R.x.value = 0;    
                R.y.value = 1;    
                R.z.value = 0;    
            } else {
                R.x = Rx;
                R.y = Ry;
                R.z = Rz;
            }
        }
    }
};
// static変数をヘッダファイルに置くと, 複数のファイルからincludeされるときにリンクエラーが起きる.
//Fp EllipticCurve::a;
//Fp EllipticCurve::b;


void add(Point& R, const Point& P, const Point& Q) {

    if (P.z.value == 0) {
        R.x = Q.x;
        R.y = Q.y;
        R.z = Q.z;
    } else if (Q.z.value == 0) {
        R.x = P.x;
        R.y = P.y;
        R.z = P.z;
    } else {
        Fp Rx, Ry, Rz;
        if (isEqual(P, Q)) { // if P == Q then Doubling
#if 1
            return EllipticCurve::dbl(R, P);
#else
            Fp u, v, v2, w, s, t;

            Fp::mulInt(s, P.x, 3);
            //mul(s, three, P.x); // 3*X
            mul(s, s, P.x); // 3*X^2
            mul(t, EllipticCurve::a, P.z); // a*Z
            mul(t, t, P.z); // a*Z^2
            add(u, s, t); // 3*X^2 + a*Z^2
        
            mul(v, P.y, P.z); // Y*Z

            mul(s, u, u); // u^2
            Fp::mulInt(t, P.x, 8);
            //mul(t, eight, P.x); // 8*X
            mul(t, t, P.y); // 8*X*Y
            mul(t, t, v); // 8*X*Y*v
            sub(w, s, t); // u^2 - 8*X*Y*v

            Fp::mulInt(Rx, v, 2);
            //mul(Rx, two, v); // 2*v
            mul(Rx, Rx, w); // Rx = 2*v*w

            Fp::mulInt(s, P.x, 4);
            //mul(s, four, P.x);
            mul(s, s, P.y);
            mul(s, s, v);
            sub(s, s, w); // 4*X*Y*v - w
            mul(s, s, u); // u(4*X*Y*v - w)

            mul(v2, v, v); // v^2

            mul(t, P.y, P.y);
            mul(t, t, v2); // (Yv)^2
            Fp::mulInt(t, t, 8);
            //mul(t, t, eight); // 8(Yv)^2

            sub(Ry, s, t); // Ry = u(4*X*Y*v - w) - 8(Yv)^2

            Fp::mulInt(Rz, v, 8);
            //mul(Rz, eight, v); // 8*v
            mul(Rz, Rz, v2); // Rz = 8*v^3
#endif
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

            Fp::mulInt(t, v2, 2);
            //mul(t, two, v2); // 2 * v^2
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

        if (Rz.value == 0) {
            R.x.value = 0;    
            R.y.value = 1;    
            R.z.value = 0;    
        } else {
            R.x = Rx;
            R.y = Ry;
            R.z = Rz;
        }
    }
}

void sub(Point& R, const Point& P, const Point& Q) {
    Point minus_Q = Q;
    sub(minus_Q.y, Fp::modulus, minus_Q.y); // y <- p - y
    add(R, P, minus_Q);
}

void mul(Point& R, const Point& P, const mpz_class& x) {
#if 0
    Fp one, zero;
    one = Fp(1);
    zero = Fp(0);

    Point k = Point(zero, one, zero);
#else
    Point k = Point(0, 1, 0);
#endif
    Point tmp_P = P;

    mpz_class n = x;
    while (n > 0) {
        if ((n % 2) == 1) {
            add(k, k, tmp_P);
        }
        //add(tmp_P, tmp_P, tmp_P);
        EllipticCurve::dbl(tmp_P, tmp_P);
        n >>= 1;
        
    }
    R.x = k.x;
    R.y = k.y;
    R.z = k.z;
}

bool isEqual(const Point& P, const Point& Q) {
    Fp s, t, u, v;

    mul(s, P.x, Q.y); // X  * Y'
    mul(t, Q.x, P.y); // X' * Y
    
    mul(u, P.x, Q.z); // X  * Z'
    mul(v, Q.x, P.z); // X' * Z
 
    return (s == t) && (u == v);
}

void print(const Point& P) {
    if (P.z.value == 0) {
        std::cout << "(0 : 1 : 0)" << std::endl;
    } else {
        Fp x_z, y_z;
        P.xy(x_z, y_z);
        std::cout << "(" << x_z.value << " : " << y_z.value << " : 1)" << std::endl;
    }
}

// Attacks against weak elliptic curve
//
void pollard_rho_f(const Point& alpha, const Point& beta, Point& x,
        mpz_class& a, mpz_class& b, const mpz_class& order) {
    Fp x_x, x_y;
    x.xy(x_x, x_y);    

    if (x_x.value % 3 == 0) {
        x = beta + x;

        b = (b + 1) % order;
    } else if (x_x.value % 3 == 1) {
        x = x * 2;
        a = (a * 2) % order;
        b = (b * 2) % order;
    } else if(x_x.value % 3 == 2) {
        x = alpha + x;
        a = (a + 1) % order;

    }
}

mpz_class pollard_rho_ECDLP(const Point& alpha, const Point& beta, 
        const EllipticCurve& curve, const mpz_class& order) {
    // beta = [d] * alpha
    mpz_class d = 0;
    mpz_class a, b, A, B;

    a = 0; b = 0;
    A = 0; B = 0;
    Point x = curve(0, 1, 0);
    Point X = curve(0, 1, 0);
    
    for(int i = 0; i < order; i++) {
        pollard_rho_f(alpha, beta, x, a, b, order);
        pollard_rho_f(alpha, beta, X, A, B, order);
        pollard_rho_f(alpha, beta, X, A, B, order);

        if (x == X && !(A == 0 && a == 0 && b == 0 && B == 0)) {
            d = (B - b) % order;
            mpz_invert(d.get_mpz_t(), d.get_mpz_t(), order.get_mpz_t());
            d = (d * (a - A)) % order;
            if (d < 0) {
                d = d + order;
            }
            return d;
        }
    }
    return d;
}


