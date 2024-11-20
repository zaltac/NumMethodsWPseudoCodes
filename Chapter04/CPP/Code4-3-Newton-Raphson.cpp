#include <cmath>
#include <iostream>
#include <iomanip>

// ==========================================================================          
// User-defined function providing f(x), which should be cast as func(x)=0.
// ========================================================================== 
double func(double x) {
    return 4.0 + x * x * (8.0 - x * x);
}

// ==========================================================================          
// User-defined function providing f'(x)=funcp(x).
// ========================================================================== 
double funcp(double x) {
    return 4.0 * x * (4.0 - x * x);
}

// ==================================================================================
// CODE4.3-NEWTON_RAPHSON.CPP. A C++ module implementing Pseudocode 4.3.                          
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
// DESCRIPTION: A C++ module to compute a root of a nonlinear equation using the                   
//   Newton-Raphson method.                                                                    
//                                                                                             
// ON ENTRY                                                                                    
//   root  :: Initial guess for the root;                                                      
//   maxit :: Maximum number of iterations permitted;                                          
//   eps   :: Convergence tolerance.                                                           
//                                                                                             
// ON EXIT                                                                                     
//   iter  :: Number of iterations realized;                                                   
//   root  :: Computed approximation for the root.                                             
//                                                                                             
// USES                                                                                        
//   ABS   :: Built-in Intrinsic function returning the absolute value of a real value;        
//                                                                                             
// ALSO REQUIRED                                                                               
//   FUNC  :: User-defined external function providing the nonlinear equation, f(x).           
//   FUNCP :: User-defined external function providing the first derivative                    
//            of the nonlinear equation, f'(x).                                                
//                                                                                             
// REVISION DATE :: 11/20/2024                                                                 
// ==================================================================================
void newton_raphson(double& root, int maxit, double eps, int& iter) {
    double small, fn, fpn, aerr, rate, x0, xn, del, del0;

    small = 1.0e-15;

    std::cout << std::endl;

    del0 = 1.0;
    x0 = root;
    iter = 0;

    do {
        fn = func(x0);
        fpn = funcp(x0);
        del = -fn / fpn;
        aerr = std::fabs(del);
        rate = aerr / (del0 * del0);
        std::cout << iter << " " << x0 << " " << fn << " " << fpn << " " << aerr << " " << rate << std::endl;
        xn = x0 + del;
        x0 = xn;
        del0 = std::fabs(del);
        iter++;
    } while ((std::fabs(fn) > eps && aerr > eps) && (iter <= maxit));

    root = xn;

    if (iter >= maxit) {
        std::cout << "** Max iteration number reached=" << std::endl;
        std::cout << "** Estimated root has NOT converged, del, f(x)= " << std::fabs(del) << " " << func(xn) << std::endl;
    }
}


// ==============================================================================
//  Main program to test module NEWTON_RAPHSON
// ==============================================================================
int main() {
    int maxit;
    double eps, root;
    int iter;

    maxit = 99;
    eps = 0.50e-6;
    root = 2.20;

    newton_raphson(root, maxit, eps, iter);

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Root is " << root << " converged after " << iter << " iterations" << std::endl;
    std::cout << "-----------------------------" << std::endl;

    return 0;
}