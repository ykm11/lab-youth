#include "curve.h"
#include "FP.h"

#include<gmpxx.h>

Fp EllipticCurve::a;
Fp EllipticCurve::b;

Fp GLV::rw;
Point GLV::base;
mpz_class GLV::lmd;
mpz_class GLV::order;


void EllipticCurve::dbl(Point &R, const Point &P) {
    if (P.z.value == 0) {
        R.x.value = 0;
        R.y.value = 1;
        R.z.value = 0;
        return;
    } 
    Fp Rx, Ry, Rz;
    Fp u, v, w, s, t;

    sqr(s, P.x); // X^2
    Fp::mulInt(s, s, 3); // 3*X^2
    sqr(t, P.z); // Z^2
    mul(t, t, a); // a*Z^2
    add(u, s, t); // 3*X^2 + a*Z^2

    mul(v, P.y, P.z); // Y*Z

    mul(w, P.x, P.y);

    sqr(s, u); // u^2
    Fp::mulInt(t, w, 8); // 8 * X * Y
    mul(t, t, v); // 8*X*Y*v
    sub(t, s, t); // w := u^2 - 8*X*Y*v

    Fp::mulInt(Rx, v, 2); // 2*v
    mul(Rx, Rx, t); // Rx = 2*v*w

    Fp::mulInt(s, w, 4);
    mul(s, s, v);
    sub(s, s, t); // 4*X*Y*v - w
    mul(s, s, u); // u(4*X*Y*v - w)


    sqr(w, v); // v^2

    sqr(t, P.y); // Y^2
    mul(t, t, w); // (Yv)^2
    Fp::mulInt(t, t, 8); // 8(Yv)^2

    sub(Ry, s, t); // Ry = u(4*X*Y*v - w) - 8(Yv)^2

    Fp::mulInt(Rz, v, 8); // 8*v
    mul(Rz, Rz, w); // Rz = 8*v^3

    if (Rz.value == 0) {
        R.x.value = 0;
        R.y.value = 1;
        R.z.value = 0;
        return;
    } 

    R.x.value = std::move(Rx.value);
    R.y.value = std::move(Ry.value);
    R.z.value = std::move(Rz.value);
}

void add(Point &R, const Point &P, const Point &Q) {

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

    mul(w, P.z, Q.z); // Z1 * Z2
    mul(Rz, v3, w); // Rz = v^3 * Z1 * Z2

    sqr(s, u); // u^2
    mul(s, s, w); // u^2 * Z1 * Z2

    Fp::mulInt(w, t, 2); // 2 * X1 * Z2 
    mul(w, w, v2); // 2 * v^2 * X1 * Z2

    sub(w, s, w); //  (u^2 * Z1 * Z2) - (2 * v^2 * X1 * Z2)
    sub(w, w, v3); // (u^2 * Z1 * Z2) - v^3 - (2 * v^2 * X1 * Z2)

    mul(Rx, v, w); // Rx = v * w

    mul(s, v2, t); // v^2 * X1 * Z2
    mul(t, v3, P.y); // v^3 * Y1
    mul(t, t, Q.z); // v^3 * Y1 * Z2
    sub(s, s, w); // v^2 * X1 * Z2 - w
    mul(s, s, u); // u(v^2 * X1 * Z2 - w)
    sub(Ry, s, t); // Ry = u(v^2 * X1 * Z2 - w) -  v^3 * Y1 * Z2

    if (Rz.value == 0) {
        R.x.value = 0;
        R.y.value = 1;
        R.z.value = 0;
        return;
    }

    R.x.value = std::move(Rx.value);
    R.y.value = std::move(Ry.value);
    R.z.value = std::move(Rz.value);
}

void sub(Point &R, const Point &P, const Point &Q) {
    Point::neg(R, Q); // R <- [-1]Q
    add(R, P, R); // R <- P + [-1]Q
}

bool isEqual(const Point &P, const Point &Q) {
    Fp s, t, u, v;

    mul(s, P.x, Q.y); // X  * Y'
    mul(t, Q.x, P.y); // X' * Y

    mul(u, P.x, Q.z); // X  * Z'
    mul(v, Q.x, P.z); // X' * Z

    return (s == t) && (u == v);
}

void dump(const Point &P) {
    if (P.z.value == 0) {
        std::cout << "(0 : 1 : 0)" << std::endl;
    } else {
        Fp x, y;
        P.xy(x, y);
        std::cout << "(" << x.value << " : " << y.value << " : 1)" << std::endl;
    }
}

void mul(Point &R, const Point &P, const mpz_class &x) { // 左向きバイナリ法
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

void r_mul(Point &R, const Point& G, const mpz_class &x) { // 右向きバイナリ法
    size_t k_bits = mpz_sizeinbase(x.get_mpz_t(), 2);

    R = G;
    for (int i = k_bits-2; i >= 0; i--) {
        EllipticCurve::dbl(R, R);

        if (mpz_tstbit(x.get_mpz_t(), i) == 1) {
            add(R, R, G);
        }
    }
}

void montgomery_mul(Point &R0, const Point &G, const mpz_class &n) { // Montgomery Ladder
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

void window_mul(Point &R, const Point &G, const mpz_class &n) { // windows method
    size_t k_bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    Point P[4];

    P[0] = Point(0, 1, 0);
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

void multipleMul(Point &R, const Point &P, const mpz_class &u, const Point &Q, const mpz_class &v) {
    size_t k_bits;
    Point prePoints[4][4]; // w = 2
    prePoints[0][0] = Point(0, 1, 0);
    prePoints[1][0] = P;
    prePoints[0][1] = Q;
    EllipticCurve::dbl(prePoints[2][0], P);
    EllipticCurve::dbl(prePoints[0][2], Q);
    add(prePoints[3][0], prePoints[2][0], P);
    add(prePoints[0][3], prePoints[0][2], Q);
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            add(prePoints[i][j], prePoints[i][0], prePoints[0][j]); // [i]P + [j]Q
        }
    }

    R = prePoints[0][0]; 
    if (mpz_sizeinbase(u.get_mpz_t(), 2) > mpz_sizeinbase(v.get_mpz_t(), 2)) {
        k_bits = mpz_sizeinbase(u.get_mpz_t(), 2);
    } else {
        k_bits = mpz_sizeinbase(v.get_mpz_t(), 2);
    }
    for (int i = k_bits-1; i > 0; i=i-2) {
        EllipticCurve::dbl(R, R);
        EllipticCurve::dbl(R, R); // R <- [4]R
        add(R, R, prePoints
                [2*mpz_tstbit(u.get_mpz_t(), i) + mpz_tstbit(u.get_mpz_t(), i-1)]
                [2*mpz_tstbit(v.get_mpz_t(), i) + mpz_tstbit(v.get_mpz_t(), i-1)]);
        // R <- R + ([]P + []Q)
    }
    if ((k_bits & 1) == 1) { // ビット数が奇数のときだけ
        EllipticCurve::dbl(R, R);
        add(R, R, prePoints
                [mpz_tstbit(u.get_mpz_t(), 0)]
                [mpz_tstbit(v.get_mpz_t(), 0)]);
    }
}

void naf_mul(Point &R, const Point &P, const mpz_class &x) {
    size_t naf_size = mpz_sizeinbase(x.get_mpz_t(), 2) + 1;
    int8_t naf[naf_size];
    memset(naf, 0, naf_size);
#if 0
    size_t w_size = 1;
    size_t tblSize = 1 << w_size;
    Point tbl[tblSize];
    tbl[0] = Point(0, 1, 0);
    tbl[1] = P;
    
    getNafArray(naf, x);
    while (naf_size >= 1 && naf[naf_size-1] == 0) {
        naf_size--;
    }

    Point Q;
    R = P;
    int8_t t;
    for (int i = naf_size-2; i >= 0; i--) {
        EllipticCurve::dbl(R, R);
        t = naf[i]; 
        if (t < 0) {
            Point::neg(Q, tbl[-t]);
        } else {
            Q = tbl[t];
        }
        add(R, R, Q);
    }
#else
    size_t w_size = 3;
    size_t tblSize = 1 << w_size;
    Point tbl[tblSize];
    tbl[0] = Point(0, 1, 0);
    tbl[1] = P;
    EllipticCurve::dbl(tbl[2], P);
    add(tbl[3], tbl[2], P);
    EllipticCurve::dbl(tbl[4], tbl[2]);
    add(tbl[5], tbl[4], P);
    
    getNafArray(naf, x);
    while (naf_size >= 1 && naf[naf_size-1] == 0) {
        naf_size--;
    }

    Point Q;
    R = P;
    int8_t t;
    int i;
    for (i = naf_size-2; i >= 2; i=i-3) {
        EllipticCurve::dbl(R, R);
        EllipticCurve::dbl(R, R);
        EllipticCurve::dbl(R, R);

        t = 4*naf[i] + 2*naf[i-1] + naf[i-2]; 
        if (t < 0) {
            Point::neg(Q, tbl[-t]);
        } else {
            Q = tbl[t];
        }
        add(R, R, Q);
    }

    if (i < 0) {
        return;
    }

    t = 0;
    if (i == 1) {
        EllipticCurve::dbl(R, R);
        t = 2*naf[1];
    }
    EllipticCurve::dbl(R, R);
    t = t + naf[0];

    if (t < 0) {
        Point::neg(Q, tbl[-t]);
    } else {
        Q = tbl[t];
    }
    add(R, R, Q);
#endif
}


void GLV::decomposing_k(mpz_class &k0, mpz_class &k1, const mpz_class &k) {
    // k = k0 + k1*lmd
    // EEAの過程で見つける.
    mpz_class t1, t0, n;
    mpz_class q, r;
    mpz_class b1, b2;

    t0 = 0; t1 = 1;
    b2 = GLV::order; n = GLV::lmd;
    mpz_sqrt(b1.get_mpz_t(), GLV::order.get_mpz_t()); // sqrt(|E|)
    while (n != 0) {
        mpz_tdiv_qr(q.get_mpz_t(), r.get_mpz_t(), b2.get_mpz_t(), n.get_mpz_t());
        k0 = t1; t1 = t0 - q*t1; t0 = k0;

        b2 = n;
        n = r;
        if (b2 < b1) {
            k0 = r; 
            k1 = -t1;
            mpz_tdiv_qr(q.get_mpz_t(), r.get_mpz_t(), b2.get_mpz_t(), n.get_mpz_t());
            t1 = q*t1 - t0;
            break;
        }
    }
    b2 = (k1*k) / (-k0*t1 + r*k1);
    b1 = -b2*t1/k1;

    t0 = b1*k0 + b2*r; // x
    r = b1*k1 + b2*t1; // y

    /*
        k0 = k - x
        k1 = -y
    */
    mpz_sub(k0.get_mpz_t(), k.get_mpz_t(), t0.get_mpz_t());
    mpz_neg(k1.get_mpz_t(), r.get_mpz_t());
}

void GLV::lambdaMul(Point &R, const Point &P) { 
    mulMod(R.x.value, GLV::rw.value, P.x.value, Fp::modulus); 
    R.y = P.y;
    R.z = P.z;
}

void GLV::mulBase(Point &R, const mpz_class &k) { 
    GLV::scalarMul(R, GLV::base, k);
}

void GLV::scalarMul(Point &R, const Point &P, const mpz_class &k) { 
    mpz_class k0, k1; 

    /*
        k = k0 + k1 * lmd
    */
    decomposing_k(k0, k1, k); 

    size_t k_bits;
    const size_t w = 2; 
    const size_t tblSize = 1 << w; // w = 2^{w_size}
    Point prePoints[tblSize][tblSize]; 
    prePoints[0][0] = Point(0, 1, 0);
    prePoints[1][0] = P;
    EllipticCurve::dbl(prePoints[2][0], P);
    for (size_t i = 2; i < tblSize-1; i++) {
        add(prePoints[i+1][0], prePoints[i][0], P);
    }
    for (size_t i = 1; i < tblSize; i++) { // Q <- [i * lmd]P
        GLV::lambdaMul(prePoints[0][i], prePoints[i][0]);
    }
    for (size_t i = 1; i < tblSize; i++) {
        for (size_t j = 1; j < tblSize; j++) {
            add(prePoints[i][j], prePoints[i][0], prePoints[0][j]); // [i]P + [j]Q
        }
    }

    R = prePoints[0][0]; 
    if (mpz_sizeinbase(k0.get_mpz_t(), 2) > mpz_sizeinbase(k1.get_mpz_t(), 2)) {
        k_bits = mpz_sizeinbase(k0.get_mpz_t(), 2);
    } else {
        k_bits = mpz_sizeinbase(k1.get_mpz_t(), 2);
    }
    for (int i = k_bits-1; i > 0; i=i-2) {
        EllipticCurve::dbl(R, R);
        EllipticCurve::dbl(R, R); // R <- [4]R
        add(R, R, prePoints
                [2*mpz_tstbit(k0.get_mpz_t(), i) + mpz_tstbit(k0.get_mpz_t(), i-1)]
                [2*mpz_tstbit(k1.get_mpz_t(), i) + mpz_tstbit(k1.get_mpz_t(), i-1)]);
        // R <- R + ([]P + []Q)
    }

    if ((k_bits & 1) == 1) { // ビット数が奇数のときだけ
        EllipticCurve::dbl(R, R);
        add(R, R, prePoints
                [mpz_tstbit(k0.get_mpz_t(), 0)]
                [mpz_tstbit(k1.get_mpz_t(), 0)]);
    }
}

