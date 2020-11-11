#include <iostream>
#include <kv/bessel.hpp>
typedef kv::interval<double> itv;
#define PRECISION 1e-13 // Set your required precision.
int nu = 0; // Set the order of the Bessel function (nu=0 or 1 is currently acceptable).

int main() {
	std::cout.precision(17);
    double eps = pow(2., -4.);
    itv x = itv(1, 1 + eps);
    itv J;
    
    std::cout << "Step 1: Preparation for bisection." << "\n";
    while (1) {
        J = kv::besselj(nu, x);
        std::cout << "x =" << x << "\n";
        std::cout << "J =" << J << "\n";
        if(J.lower()*J.upper() <= 0.) {
            std::cout << "Bisection staring from here." << "\n";
            break;
        }
        x += eps/2;
    }
    std::cout << "============================" << "\n";
    itv a, b, s, Ja, Jb, Js, JaJb, JaJs, JbJs;
    std::cout << "Step 2: Bisection method" << "\n";
    a = x.lower();
    b = x.upper();
    
    //Check sup(JaJb)<0
    Ja = kv::besselj(nu, a);
    Jb = kv::besselj(nu, b);
    JaJb = Ja * Jb;
    if (JaJb.upper() >= 0) {
        std::cout << "Error: Verification failed." <<"\n";
        return 0;
    }
    
    //Bisection method
    while (1) {
        std::cout << "a=" << a << " b=" << b << "\n";
        std::cout << "Ja=" << Ja << " Jb=" << Jb << "\n";
        if (mag(Ja - Jb) < PRECISION) {
            std::cout <<"The first zero of J_"<< nu << " is in " << itv::hull(a, b) << "\n";
            break;
        }
        s = (a + b) / 2.;
        Js = kv::besselj(nu, s);
        JaJs = Ja * Js;
        JbJs = Jb * Js;
        if (JaJs.upper() < 0) {
            b = s;
            Jb = Js;
        }
        else if(JbJs.upper() < 0) {
            a = s;
            Ja = Js;
        }
        else{
            std::cout <<"Precision can be no more refined." << "\n";
            std::cout <<"The first zero of J_"<< nu << " is in " << itv::hull(a, b) << "\n";
            break;
        }
    
    }
}
