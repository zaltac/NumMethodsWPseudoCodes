#include <cmath>
#include <iostream>
#include <iomanip>

double func(double x) {
// ==========================================================================          
// User-defined function providing f(x), which should be cast as Func(x)=0.
// ==========================================================================  
    return x * x + 0.025 * x - 4.0;
}

// ==================================================================================
// CODE4.1-BISECTION.CPP. A C++ module implementing Pseudocode 4.1.                                   
//
// NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
// First Edition. (c) By Zekeriya ALTAÇ (2024).
// ISBN: 978-1-032-75474-1 (hbk)
// ISBN: 978-1-032-75642-4 (pbk)
// ISBN: 978-1-003-47494-4 (ebk)
// 
// DOI : 10.1201/9781003474944
// C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
// 
// This free software is complimented by the author to accompany the textbook.
// E-mail: altacz@gmail.com.
// 
// DESCRIPTION: A C++ module to find a root of a nonlinear equation in [a,b]                       
//   using the Bisection method.                                                               
//                                                                                             
// ON ENTRY                                                                                    
//  [a,b] :: Initial search interval (it must bracket one root);                               
//  maxit :: Maximum number of iterations permitted;                                           
//  eps   :: Convergence tolerance.                                                            
//                                                                                             
// ON EXIT                                                                                     
//  halves:: Number of halves realized;                                                        
//  root  :: Computed approximation for the root.                                              
//                                                                                             
// USES                                                                                        
//  ABS   :: Built-in Intrinsic function returning the absolute value of a real value;         
//                                                                                             
// ALSO REQUIRED                                                                               
//  FUNC  :: User-defined external function providing the nonlinear equation.                  
//                                                                                             
// REVISION DATE :: 11/20/2024                                                                 
// ==================================================================================
void Bisection(double a, double b, int maxit, double eps, double& root, int& halves) {
    int p = 0;
    double fa = func(a);
    double fb = func(b);
    double interval = b - a;

    std::cout << "  p     a           b           f(a)        f(b)        xm          f(xm)       interval" << std::endl;
    std::cout << "----------------------------------------------------------------------------------------------" << std::endl;

    while (true) {
        p++;
        double xm = 0.5 * (a + b);
        double fm = func(xm);
        std::cout << std::setw(3) << p << " " << std::scientific << std::setprecision(7) << a << " " << b << " " << fa << " " << fb << " " << xm << " " << fm << " " << 0.5 * interval << std::endl;

        if (fa * fm > 0.0) {
            a = xm;
            fa = fm;
        } else {
            b = xm;
            fb = fm;
        }

        interval *= 0.5;

        if ((std::fabs(fm) < eps && interval < eps) || p == maxit)
            break;
    }

    root = 0.5 * (a + b);
    halves = p;

    if (p == maxit)
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
                  << "! Max iteration number reached=" << maxit << " !" << std::endl
                  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
}

// ==================================================================================
// Main program to test BISECTION.CPP
// ==================================================================================
int main() {
    int Halves, maxit;
    double a, b, eps, root, fa, fb;

    maxit = 99;
    eps = 0.00005;
    a = 0.0;
    b = 4.0;

    fa = func(a);
    fb = func(b);

    if (fa * fb > 0) {
        std::cout << "No root in interval (a,b). Change the interval." << std::endl;
        return 1;
    }

    Bisection(a, b, maxit, eps, root, Halves);

    std::cout << "=======================================" << std::endl
              << "! Root is " << std::scientific << std::setprecision(8) << root << " after " << Halves << " bisections !" << std::endl
              << "=======================================" << std::endl
              << std::endl
              << "*** Root (with linear interpolation)   = " << std::scientific << std::setprecision(9) << (a * fb - b * fa) / (fb - fa) << std::endl
              << "*** Clossness to the root, Abs[f(root)]= " << std::scientific << std::setprecision(4) << std::fabs(func(root)) << std::endl;

    return 0;
}