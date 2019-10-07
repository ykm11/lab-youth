#include "curve.h"
#include "FP.h"

#include<gmpxx.h>

Fp EllipticCurve::a;
Fp EllipticCurve::b;

void EllipticCurve::dbl(Point& R, const Point& P) {
    if (P.z.value == 0) {
        R.x.value = 0;
        R.y.value = 1;
        R.z.value = 0;
    } else {
        Fp Rx, Ry, Rz;
        Fp u, v, v2, w, s, t;

        sqr(s, P.x); // X^2
        Fp::mulInt(s, s, 3); // 3*X^2
        sqr(t, P.z); // Z^2
        mul(t, t, a); // a*Z^2
        add(u, s, t); // 3*X^2 + a*Z^2

        mul(v, P.y, P.z); // Y*Z

        sqr(s, u); // u^2
        Fp::mulInt(t, P.x, 8); // 8*X
        mul(t, t, P.y); // 8*X*Y
        mul(t, t, v); // 8*X*Y*v
        sub(w, s, t); // u^2 - 8*X*Y*v

        Fp::mulInt(Rx, v, 2); // 2*v
        mul(Rx, Rx, w); // Rx = 2*v*w

        Fp::mulInt(s, P.x, 4);
        mul(s, s, P.y);
        mul(s, s, v);
        sub(s, s, w); // 4*X*Y*v - w
        mul(s, s, u); // u(4*X*Y*v - w)

        sqr(v2, v); // v^2

        sqr(t, P.y); // Y^2
        mul(t, t, v2); // (Yv)^2
        Fp::mulInt(t, t, 8); // 8(Yv)^2

        sub(Ry, s, t); // Ry = u(4*X*Y*v - w) - 8(Yv)^2

        Fp::mulInt(Rz, v, 8); // 8*v
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

void add(Point& R, const Point& P, const Point& Q) {

    if (P.z.value == 0) {
        R.x = Q.x;
        R.y = Q.y;
        R.z = Q.z;
        return;
    }
    if (Q.z.value == 0) {
        R.x = P.x;
        R.y = P.y;
        R.z = P.z;
        return;
    }
    Fp u, v, s, t;

    mul(s, Q.y, P.z); // Y2 * Z1
    mul(t, P.y, Q.z); // Y1 * Z2
    sub(u, s, t); // Y2 * Z1 - Y1 * Z2

    mul(s, Q.x, P.z); // X2 * Z1
    mul(t, P.x, Q.z); // X1 * Z2
    sub(v, s, t); // X2 * Z1 - X1 * Z2

    if (u.value == 0 && v.value == 0) {
        EllipticCurve::dbl(R, P);
        return;
    }
    // otherwise, Adding
    Fp w, v2, v3;
    Fp Rx, Ry, Rz;

    sqr(v2, v); // v^2
    mul(v3, v2, v); // v^3

    Fp::mulInt(t, v2, 2); // 2 * v^2
    mul(t, t, P.x); // 2 * v^2 * X1
    mul(t, t, Q.z); //  2 * v^2 * X1*Z2

    sqr(s, u); // u^2
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

    if (Rz.value == 0) {
        R.x.value = 0;
        R.y.value = 1;
        R.z.value = 0;
        return;
    }
#if 0
    R.x = Rx;
    R.y = Ry;
    R.z = Rz;
#else
    R.x.value = std::move(Rx.value);
    R.y.value = std::move(Ry.value);
    R.z.value = std::move(Rz.value);
#endif
}

void sub(Point& R, const Point& P, const Point& Q) {
    Point minus_Q = Q;
    minus_Q.y.value = (Fp::modulus - minus_Q.y.value); // y <- p-y
    add(R, P, minus_Q); // R <- P + [-1]Q
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

void mul(Point& R, const Point& P, const mpz_class& x) { // 左向きバイナリ法
    R.x.value = 0;
    R.y.value = 1;
    R.z.value = 0;

    Point tmp_P = P;

    mpz_class n = x;
    while (n > 0) {
        if ((n & 1) == 1) {
            add(R, R, tmp_P);
        }
        EllipticCurve::dbl(tmp_P, tmp_P);
        n >>= 1;
    }
}

void r_mul(Point &R, const Point& G, const mpz_class x) { // 右向きバイナリ法
    size_t k_bits = mpz_sizeinbase(x.get_mpz_t(), 2);

    R = G;
    for (int i = k_bits-2; i >= 0; i--) {
        EllipticCurve::dbl(R, R);

        //if (((n >> i) & 1) == 1) {
        if (mpz_tstbit(x.get_mpz_t(), i) == 1) {
            add(R, R, G);
        }
    }
}

void montgomery_mul(Point &R0, const Point& G, const mpz_class n) { // Montgomery Ladder
    size_t k_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    Point R1 = G;
    R0 = Point(0, 1, 0);

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

void window_mul(Point &R, const Point& G, const mpz_class n) { // windows method
    size_t k_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    Point P[4];

    P[0] = Point(0, 1, 0);
    P[1] = G;
    EllipticCurve::dbl(P[2], G);
    add(P[3], P[2], G);

    R = P[2*mpz_tstbit(n.get_mpz_t(), k_bits-1) + mpz_tstbit(n.get_mpz_t(), k_bits-2)];
    for (int i = k_bits-3; i > 0; i=i-2) {
        EllipticCurve::dbl(R, R);
        EllipticCurve::dbl(R, R); // R <- 4R

        add(R, R, P[2*mpz_tstbit(n.get_mpz_t(), i) + mpz_tstbit(n.get_mpz_t(), i-1)]);
    } 
    if ((k_bits & 1) == 1) { // nのビット数が奇数のときだけ
        EllipticCurve::dbl(R, R);
        add(R, R, P[mpz_tstbit(n.get_mpz_t(), 0)]);
    }
}
