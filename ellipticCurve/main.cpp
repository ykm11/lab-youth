#include<iostream>
#include "curve.h"

int main() {
    Fp::setModulo(19);

    EllipticCurve EC = EllipticCurve(2, 3);

    //std::cout << "[+] Elliptic Curve" << std::endl;
    //std::cout << EC.a.value << std::endl;

    Point P = EC.point(1, 5);
    Point Q = EC.point(3, 6);
    std::cout << P.x.value << " " << P.y.value << " " << P.z.value << std::endl;
    std::cout << Q.x.value << " " << Q.y.value << " " << Q.z.value << std::endl;

    Point R = P + Q;
    std::cout << R.x.value << " " << R.y.value << " " << R.z.value << std::endl;

}
