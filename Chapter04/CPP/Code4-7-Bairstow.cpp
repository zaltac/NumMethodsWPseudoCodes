#include <iostream>
#include <cmath>
#include <iomanip>


// ==================================================================================
// CODE4.7-BAIRSTOW.CPP. A C++ module implementing Pseudocode 4.7.                                    
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
// DESCRIPTION: A C++ module to find all real and/or imaginary roots of a polynomial           
//   of the n'th degree using the BAIRSTOW's method.                                         
//                                                                                             
// ON ENTRY                                                                                    
//   n    :: Degree of the polynomial;                                                         
//  p0,q0 :: Initial guesses for a quadratic equation; i.e., for p and q;                      
//    a   :: Array of length (n+1) containing the coefficients of polynomial defined as        
//                a0 x^n + a1 x^(n-1) + ... + an = 0                                           
//   eps  :: Convergence tolerance;                                                            
//  maxit :: Maximum number of iterations permitted;                                           
//  iprnt :: printing key, =0 do not print intermediate results, <> 0 print intermediates.     
//                                                                                             
// ON EXIT                                                                                     
//   xre  :: Array of length n containing real parts of the roots;                             
//   xim  :: Array of length n containing imaginary parts of the roots.                        
//                                                                                             
// OTHER VARIABLES                                                                             
//    b   :: Array of length [n] containing coefficients of quotient polynomial (0<=k<=n-2);   
//    c   :: Array of length [n] containing coefficients of partial derivatives.               
//                                                                                             
// USES                                                                                        
//   fabs :: Built-in Intrinsic function returning the absolute value of a real value;         
//   QUADRATIC :: A module that solves a quadratic equation of the form x2 + p x + q = 0. (see CODE1-3)    
//                                                                                             
// REVISION DATE :: 04/29/2024                                                                 
// ==================================================================================
void BAIRSTOW(int n, double p0, double q0, double a[], double eps, int maxit, int iprnt, double xre[], double xim[]) {
    double b[20], c[20], xr[2], xi[2];
    double p, delp, q, delq, delM, cbar, del, del1, del2;
    int i, k, m, kount;

    for (k = n; k >= 0; k--) {
        a[k] /= a[0]; // Normalize a's by a(0)
    }
    m = n; // Save n for later use
    kount = 0;

    while (n > 1) {
        k = 0;
        p = p0;
        q = q0; // Initialize
        delM = 1.0;

        while (delM > eps && k <= maxit) { // inner loop
            k++;
            b[0] = 1.0;
            c[0] = 1.0;
            b[1] = a[1] - p;
            c[1] = b[1] - p;

            for (i = 2; i <= n; i++) {
                b[i] = a[i] - p * b[i - 1] - q * b[i - 2];
                c[i] = b[i] - p * c[i - 1] - q * c[i - 2];
            }

            cbar = c[n - 1] - b[n - 1];
            del = c[n - 2] * c[n - 2] - cbar * c[n - 3];
            del1 = b[n - 1] * c[n - 2] - b[n] * c[n - 3];
            del2 = b[n] * c[n - 2] - b[n - 1] * cbar;
            delp = del1 / del;
            delq = del2 / del;
            p += delp;
            q += delq;
            delM = std::fabs(delp) + std::fabs(delq);

            if (iprnt == 1) {
                std::cout << "Iter= " << k << "   delM= " << delM << "   p= " << p << "   q= " << q << std::endl;
            } else if (iprnt == 2) {
                std::cout << "\niter=" << k << std::endl;
                std::cout << "---------" << std::endl;
                std::cout << " dp =" << delp << "   dq =" << delq << "   delM=" << delM << std::endl;
                std::cout << "  p =" << p << "    q =" << q << std::endl;
                std::cout << "-------------------------------------" << std::endl;
                std::cout << " k     a(k)       b(k)       c(k)" << std::endl;
                for (i = 0; i <= n; i++) {
                    std::cout << " " << i << "  " << a[i] << "  " << b[i] << "  " << c[i] << std::endl;
                }
                std::cout << "-------------------------------------" << std::endl;
            }
        } // end of inner loop

        if (k - 1 == maxit) {
            std::cout << "Quadratic factor did not converge after " << k - 1 << " iterations" << std::endl;
            std::cout << "Recent values of p, q, delM are " << p << ", " << q << ", " << delM << std::endl;
            std::cout << "Corresponding roots may be questionable ..." << std::endl;
        }

        QUADRATIC(p, q, xr, xi);
        xre[kount] = xr[0];
        xim[kount] = xi[0];
        kount++;
        xre[kount] = xr[1];
        xim[kount] = xi[1];
        kount++;

        std::cout << std::endl;
        std::cout << "======== FOUND A QUADRATIC FACTOR ========" << std::endl;
        std::cout << "   x*x + (" << p << ")*x + (" << q << ")" << std::endl;
        std::cout << "==========================================" << std::endl << std::endl;

        n -= 2;
        for (i = 0; i <= n; i++) {
            a[i] = b[i];
        }

        if (n == 1) {
            xre[kount] = -a[1];
            xim[kount] = 0.0;
            kount++;
            std::cout << "======== FOUND A LINEAR FACTOR ========" << std::endl;
            std::cout << "      x  + (" << a[1] << ")" << std::endl;
            std::cout << "======================================" << std::endl << std::endl;
        }
    }

    n = m;
}


// ==============================================================================
//  Main program to test BAIRSTOW module
// ==============================================================================
int main() {
    double a[20], xre[19], xim[19];
    int n = 5, iprnt = 2, maxit = 99;
    double p0 = 0.0, q0 = 0.0, eps = 0.5e-4;

    a[0] = 1.0;
    a[1] = -5.0;
    a[2] = -15.0;
    a[3] = 85.0;
    a[4] = -26.0;
    a[5] = -120.0;

    BAIRSTOW(n, p0, q0, a, eps, maxit, iprnt, xre, xim);

    std::cout << "    ======== All the Roots are =========" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "Root(" << i + 1 << ") = " << xre[i] << " + (" << xim[i] << ") i" << std::endl;
    }
    std::cout << "    ===================================" << std::endl;

    return 0;
}