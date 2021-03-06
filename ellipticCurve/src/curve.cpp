#include "curve.h"
#include "FP.h"

#include<gmpxx.h>

Fp EllipticCurve::a;
Fp EllipticCurve::b;

void EllipticCurve::dbl(Point &R, const Point &P) {
    if (zeroCmp(P.z)) {
        setInfPoint(R);
        return;
    } 
    Fp Rx, Ry, Rz;
    Fp u, v, w, s, t;

#ifndef YKM_ECC_USE_MPN
    s.value = a.value + 3;
    if (s.value == Fp::modulus) { // a == -3 ? 
#else
    add(s, a, 3);
    if (zeroCmp(s)) { // a == -3 ? 
#endif
        sub(s, P.x, P.z); // X - Z
        add(t, P.x, P.z); // X + Z
        mul(u, s, t); // (X - Z) * (X + Z)
        Fp::mulInt(u, u, 3); // 3 * (X - Z) * (X + Z)
    } else {
        sqr(s, P.x); // X^2
        Fp::mulInt(s, s, 3); // 3*X^2
        sqr(t, P.z); // Z^2
        mul(t, t, a); // a*Z^2
        add(u, s, t); // 3*X^2 + a*Z^2
    }

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

    if (zeroCmp(Rz)) {
        setInfPoint(R);
        return;
    } 

    setPoint(R, Rx, Ry, Rz);
}

void add(Point &R, const Point &P, const Point &Q) {

    if (zeroCmp(P.z)) {
        setPoint(R, Q.x, Q.y, Q.z);
        return;
    } else if (zeroCmp(Q.z)) {
        setPoint(R, P.x, P.y, P.z);
        return;
    }
    Fp u, v, s, t;

    mul(s, Q.y, P.z); // Y2 * Z1
    mul(t, P.y, Q.z); // Y1 * Z2
    sub(u, s, t); // Y2 * Z1 - Y1 * Z2

    mul(s, Q.x, P.z); // X2 * Z1
    mul(t, P.x, Q.z); // X1 * Z2
    sub(v, s, t); // X2 * Z1 - X1 * Z2

    if (zeroCmp(u) && zeroCmp(v)) {
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

    if (zeroCmp(Rz)) {
        setInfPoint(R);
        return;
    }

    setPoint(R, Rx, Ry, Rz);
}


void EllipticCurve::dbl(jPoint &R, const jPoint &P) {
    if (zeroCmp(P.y)) {
        setInfPoint(R);
        return;
    }

    Fp u, v, s, t;
    Fp Rx, Ry, Rz;
    sqr(u, P.y); // Y^{2}
    mul(s, u, P.x); // X * Y^{2}
    Fp::mulInt(s, s, 4); // 4 * X * Y^{2}
    
    sqr(t, P.z); // Z^{2}

#ifndef YKM_ECC_USE_MPN
    v.value = a.value + 3;
    if (v.value == Fp::modulus) {
#else
    add(v, a, 3);
    if (zeroCmp(v)) { // a == -3 ? 
#endif
        add(Rx, P.x, t); // (X + Z^2)
        sub(v, P.x, t); // (X - Z^2)
        mul(v, v, Rx); // (X + Z^2) * (X - Z^2)
        Fp::mulInt(v, v, 3); // 3 * (X + Z^2) * (X - Z^2)
    } else {
        sqr(t, t); // Z^{4}
        mul(t, t, a); // a * Z^{4}
        sqr(v, P.x); // X^{2}
        Fp::mulInt(v, v, 3); // 3 * X^{2}
        add(v, v, t); // 3 * X^{2} + a * Z^{4}
    }

    sqr(Rx, v); // (3 * X^{2} + a * Z^{4})^2
    sub(Rx, Rx, s);
    sub(Rx, Rx, s);

    sub(Ry, s, Rx); // 
    mul(Ry, Ry, v);
    sqr(u, u); // Y^{4}
    Fp::mulInt(u, u, 8); // 8 * Y^{4}
    sub(Ry, Ry, u); // 

    mul(Rz, P.y, P.z);
    Fp::mulInt(Rz, Rz, 2);

    if (zeroCmp(Rz)) {
        setInfPoint(R);
        return;
    } 

    setPoint(R, Rx, Ry, Rz);
}

void add(jPoint &R, const jPoint &P, const jPoint &Q) {
    if (zeroCmp(P.z)) {
        setPoint(R, Q.x, Q.y, Q.z);
        return;
    } else if (zeroCmp(Q.z)) {
        setPoint(R, P.x, P.y, P.z);
        return;
    }
    Fp u, v, s, t;

    sqr(s, Q.z); // Z2^{2}
    sqr(t, P.z); // Z1^{2}
    
    mul(u, P.x, s); // X1 * Z2^{2}
    mul(v, Q.x, t); // X2 * Z1^{2}

    mul(s, s, Q.z); // Z2^{3}
    mul(t, t, P.z); // Z1^{3}

    mul(s, s, P.y); // Y1 * Z2^{3}
    mul(t, t, Q.y); // Y2 * Z1^{3}

    sub(v, v, u); // X2 * Z1^{2} - X1 * Z2^{2}
    sub(t, t, s); // Y2 * Z1^{3} - Y1 * Z2^{3}

    if (zeroCmp(v)) {
        if (zeroCmp(t)) {
            EllipticCurve::dbl(R, P);
            return;
        }
        setInfPoint(R);
        return;
    }
    Fp w;
    Fp Rx, Ry, Rz;
    mul(Rz, P.z, Q.z); 
    mul(Rz, Rz, v); // Z3 =  (X2 * Z1^{2} - X1 * Z2^{2}) * Z1 * Z2

    sqr(w, v);
    mul(u, u, w); // u * (X2 * Z1^{2} - X1 * Z2^{2})^{2}
    Fp::mulInt(Ry, u, 2);
    mul(v, w, v); // (X2 * Z1^{2} - X1 * Z2^{2})^{3}

    sqr(Rx, t); 
    sub(Rx, Rx, v); // (Y2 * Z1^{3} - Y1 * Z2^{3})^{2} - (X2 * Z1^{2} - X1 * Z2^{2})^{3}
    sub(Rx, Rx, Ry);

    sub(Ry, u, Rx);
    mul(Ry, Ry, t);
    mul(v, v, s);
    sub(Ry, Ry, v);

    if (zeroCmp(Rz)) {
        setInfPoint(R);
        return;
    } 

    setPoint(R, Rx, Ry, Rz);
}

bool isEqual(const Point &P, const Point &Q) {
    Fp s, t, u, v;

    mul(s, P.x, Q.y); // X  * Y'
    mul(t, Q.x, P.y); // X' * Y

    mul(u, P.x, Q.z); // X  * Z'
    mul(v, Q.x, P.z); // X' * Z

    return isEq(s, t) && isEq(u, v);
}

bool isEqual(const jPoint &P, const jPoint &Q) {
    Fp s, t, u, v;
    sqr(s, Q.z); // Z2^{2}
    sqr(t, P.z); // Z1^{2}
    
    mul(u, P.x, s); // X1 * Z2^{2}
    mul(v, Q.x, t); // X2 * Z1^{2}

    mul(s, s, Q.z); // Z2^{3}
    mul(t, t, P.z); // Z1^{3}

    mul(s, s, P.y); // Y1 * Z2^{3}
    mul(t, t, Q.y); // Y2 * Z1^{3}

    return isEq(s, t) && isEq(u, v);
}

void multipleMul(Point &R, const Point &P, const mpz_class &u, const Point &Q, const mpz_class &v) {
    size_t k_bits;
    Point prePoints[4][4]; // w = 2
    setInfPoint(prePoints[0][0]);
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


Fp GLV::rw;
Point GLV::base;
mpz_class GLV::lmd;
mpz_class GLV::order;

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
    mul(R.x, GLV::rw, P.x);
    copy(R.y, P.y);
    copy(R.z, P.z);
}

void GLV::mulBase(Point &R, const mpz_class &k) { 
    GLV::scalarMul(R, GLV::base, k);
}

void GLV::scalarMul(Point &R, const Point &P, const mpz_class &k) { 
    const size_t w_size = 4; 
    const size_t tblSize = 1 << w_size; // w = 2^{w_size}

    mpz_class k0, k1; 
    Point tbl0[tblSize];
    Point tbl1[tblSize];

    /*
        k = k0 + k1 * lmd
    */
    decomposing_k(k0, k1, k); 

    size_t naf0_size = mpz_sizeinbase(k0.get_mpz_t(), 2) + 1;
    size_t naf1_size = mpz_sizeinbase(k1.get_mpz_t(), 2) + 1;
    int8_t naf0[naf0_size];
    int8_t naf1[naf1_size];
    memset(naf0, 0, naf0_size);
    memset(naf1, 0, naf1_size);

    tbl0[0] = Point(0, 1, 0);
    tbl0[1] = P;
    tbl1[0] = tbl0[0];
    GLV::lambdaMul(tbl1[1], P);

    for (size_t k = 2; k < 11; k++) {
        add(tbl0[k], tbl0[k-1], P);
        GLV::lambdaMul(tbl1[k], tbl0[k]);
    }

    getNafArray(naf0, k0);
    getNafArray(naf1, k1);

    Point Q;
    int8_t t0, t1;
    R = tbl0[0];
    t0 = t1 = 0;
    while (naf0_size > naf1_size) {
        t0 <<= 1;
        t0 += naf0[naf0_size-1];
        naf0_size--;
    }

    while (naf1_size > naf0_size) {
        t1 <<= 1;
        t1 += naf1[naf1_size-1];
        naf1_size--;
    }

    if (t0 < 0) {
        Point::neg(Q, tbl0[-t0]);
    } else {
        Q = tbl0[t0];
    }
    add(R, R, Q);

    if (t1 < 0) {
        Point::neg(Q, tbl1[-t1]);
    } else {
        Q = tbl1[t1];
    }
    add(R, R, Q);

    if (naf0_size == 0 && naf1_size == 0) {
        return;
    }

    t0 = t1 = 0;
    while (naf1_size % w_size != 0) {
        EllipticCurve::dbl(R, R);
        t0 <<= 1;
        t0 += naf0[naf0_size - 1];
        naf0_size--;
        t1 <<= 1;
        t1 += naf1[naf1_size - 1];        
        naf1_size--;
    }

    if (t0 < 0) {
        Point::neg(Q, tbl0[-t0]);
    } else {
        Q = tbl0[t0];
    }
    add(R, R, Q);

    if (t1 < 0) {
        Point::neg(Q, tbl1[-t1]);
    } else {
        Q = tbl1[t1];
    }
    add(R, R, Q);

    for (int i = naf1_size-1; i >= int(w_size-1); i=i-w_size) {
        t0 = 0;
        t1 = 0;
        for (size_t j = 0; j < w_size; j++) {
            EllipticCurve::dbl(R, R);
            t0 <<= 1;
            t0 += naf0[i - j];
            t1 <<= 1;
            t1 += naf1[i - j];
        }
        
        if (t0 < 0) {
            Point::neg(Q, tbl0[-t0]);
        } else {
            Q = tbl0[t0];
        }
        add(R, R, Q);

        if (t1 < 0) {
            Point::neg(Q, tbl1[-t1]);
        } else {
            Q = tbl1[t1];
        }
        add(R, R, Q);
    }

}

Fp TwistedEdwardCurve::a;
Fp TwistedEdwardCurve::d;
Fp TwistedEdwardCurve::I_;
Point TwistedEdwardCurve::base_;

void TwistedEdwardCurve::Pneg(Point &R, const Point &P) {
#if 1
    Fp::neg(R.x, P.x);
    R.y = P.y;
#else
    R.x = P.x;
    Fp::neg(R.y, P.y);
#endif
    R.z = P.z;
}

void TwistedEdwardCurve::Padd(Point &R, const Point &P, const Point &Q) {
    Fp k, l, s, t, u, v, w;
    Fp Rx, Ry, Rz;

    if (zeroCmp(P.x)) {
        R.x = Q.x;
        R.y = Q.y;
        R.z = Q.z;
        return;
    }else if (zeroCmp(Q.x)) {
        R.x = P.x;
        R.y = P.y;
        R.z = P.z;
        return;
    }

    // How to switch to DBL? 

    mul(s, P.z, Q.z); // A := Z1 * Z2
    sqr(t, s); // B := A^2 = (Z1 * Z2) ^ 2
    mul(u, P.x, Q.x); // C := X1 * X2
    mul(v, P.y, Q.y); // D := Y1 * Y2
    mul(w, d, u);
    mul(w, w, v); // E := d * (X1 * X2) * (Y1 * Y2)

    sub(k, t, w); // F := B - E
    add(l, t, w); // G := B + E

    add(Ry, P.x, P.y);
    add(Rz, Q.x, Q.y);
    mul(Rx, Ry, Rz);
    add(Rz, u, v); // (C + D)
    sub(Rx, Rx, Rz); // - C - D = - (C + D)
    mul(Rx, Rx, k);
    mul(Rx, Rx, s);

    mul(Ry, s, l);
    mul(Ry, Ry, Rz);
    
    mul(Rz, k, l); 

    R.x = Rx;
    R.y = Ry;
    R.z = Rz;
}

void TwistedEdwardCurve::Pdbl(Point &R, const Point &P) {
    Fp k, l, s, t, u, v;
    Fp Rx, Ry, Rz;

    if (zeroCmp(P.x)) {
        R.x = P.x;
        R.y = P.y;
        R.z = P.z;
        return;
    }
    // 最適化（しなさい）

    add(k, P.x, P.y); // (X + Y)
    sqr(k, k); // B := (X + Y) ^ 2
    sqr(l, P.x); // C := X ^ 2
    sqr(s, P.y); // D := Y ^ 2
    sub(t, s, l); // F := D - C = Y ^ 2 - X ^ 2
    sqr(u, P.z); // H := Z ^ 2
    Fp::mulInt(v, u, 2);
    sub(v, t, v); // J := F - 2H

    add(u, l, s); // C + D
    sub(Rx, k, u); // B - (C + D)
    mul(Rx, Rx, v); // Rx

    mul(Ry, u, t); // -Ry

    mul(Rz, t, v); // Rz

    R.x = Rx;
    Fp::neg(R.y, Ry);
    R.z = Rz;
}

void TwistedEdwardCurve::baseMult(Point &R, const mpz_class &x) { 
    scalarMul(R, base_, x);
}

void TwistedEdwardCurve::scalarMul(Point &R, const Point &P, const mpz_class &x) { 
    size_t naf_size = mpz_sizeinbase(x.get_mpz_t(), 2) + 1;
    int8_t naf[naf_size];
    memset(naf, 0, naf_size);

    const size_t w_size = 5;
    const size_t tblSize = 1 << w_size;
    Point tbl[tblSize];
    setInfPoint(tbl[0]);
    tbl[1] = P;

    for (size_t k = 2; k < 21; k=k+2) {
        TwistedEdwardCurve::Pdbl(tbl[k], tbl[k/2]);
        TwistedEdwardCurve::Padd(tbl[k+1], tbl[k], P);
    }
 
    getNafArray(naf, x);
    while (naf_size > 0) {
        if (naf[naf_size - 1]) break;
        naf_size--;
    }

    Point Q;
    int8_t t = 1;
    while((naf_size-1) % w_size != 0) {
        t <<= 1;
        t += naf[naf_size-2];
        naf_size--;
    }
    if (t < 0) {
        TwistedEdwardCurve::Pneg(R, tbl[-t]);
    } else {
        R = tbl[t];
    }

    for (int i = naf_size-2; i >= int(w_size-1); i=i-w_size) {
        t = 0;
        for (size_t j = 0; j < w_size; j++) {
            TwistedEdwardCurve::Pdbl(R, R);
            t <<= 1;
            t += naf[i-j];
        }

        if (t < 0) {
            TwistedEdwardCurve::Pneg(Q, tbl[-t]);
        } else {
            Q = tbl[t];
        }
        TwistedEdwardCurve::Padd(R, R, Q);
    }
}

void TwistedEdwardCurve::encodePoint(uint8_t buf[32], const Point &P) {
/*
def encodepoint(P):
  x = P[0]
  y = P[1]
  bits = [(y >> i) & 1 for i in range(b - 1)] + [x & 1]
  return ''.join([chr(sum([bits[i * 8 + j] << j for j in range(8)])) for i in range(b/8)])
*/
    Fp x, y;
    P.xy(x, y);
#ifdef YKM_ECC_USE_MPN
    memcpy(buf, y.value, 32);
    buf[31] &= 0x7f;
    buf[31] |= ((uint8_t)(x.value[0] & 1) << 7);
#else
#endif
}

void TwistedEdwardCurve::xrecover(Fp &x, const Fp &y) {
/*
def xrecover(y):
  xx = (y*y-1) * inv(d*y*y+1)
  x = expmod(xx,(q+3)/8,q)
  if (x*x - xx) % q != 0: x = (x*I) % q
  if x % 2 != 0: x = q-x
  return x
*/
#ifdef YKM_ECC_USE_MPN
	Fp xx, yy;
    Fp u, v;

    mp_limb_t tp[80*Fp::size_];
    Fp e(mpz_class("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe", 16));

	sqr(yy, y); // yy
    sub_1(xx.value, yy.value, Fp::size_, 1); // yy - 1
    mul(u, d, yy); // d * yy
    add_1(u.value, u.value, Fp::size_, 1); // d * yy + 1
    invmod(v, u); // inv(d * yy + 1)
    mul(xx, xx, v);
    powMod(x.value, xx.value, e.value, Fp::modulus, tp, Fp::size_);

    // if (x*x - xx) % q != 0: x = (x*I) % q
    sqr(v, x);
    if (!isEq(v, xx)) {
        mul(x, x, I_);
    }

    // if x % 2 != 0: x = q-x
    if (x.value[0] % 2 == 1) {
        sub(x, Fp::modulus, x);
    }
#else
#endif
}


void TwistedEdwardCurve::genPublicKey(uint8_t *pk, const uint8_t *sk, size_t len) {
/*
def publickey(sk):
  h = H(sk)
  a = 2**(b-2) + sum(2**i * bit(h,i) for i in range(3,b-2))
  A = scalarmult(B,a)
  return encodepoint(A)
*/
	uint8_t hash[64];
	//sha512(hash, sk, len); // len(secret key) == 32 ?
    mp_limb_t h[4];
    memcpy(h, hash, 32);

    h[0] &= 0xFFFFFFFFFFFFFFF8;
    h[3] &= 0x7FFFFFFFFFFFFFFF;
    h[3] |= 0x4000000000000000;

    mpz_t a;
    a->_mp_d = h;
    a->_mp_size = 4;
    a->_mp_alloc = 4;

    Point R;
    mpz_class aa(a);
    baseMult(R, aa);
    encodePoint(pk, R);
}


