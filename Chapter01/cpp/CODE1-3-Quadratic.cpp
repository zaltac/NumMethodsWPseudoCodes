#include <cmath>
#include <iostream>

void quadratic_eq(double p, double q, double& re1, double& im1, double& re2, double& im2) {
// ==================================================================================
// CODE1.3-Quadratic_Eq.cpp. A fortran module for implementing Pseudocode 1.3.          
//
// NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
// First Edition. (c) By Zekeriya ALTAC (2024).
// ISBN: 978-1-032-75474-1 (hbk)
// ISBN: 978-1-032-75642-4 (pbk)
// ISBN: 978-1-003-47494-4 (ebk)
// 
// DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
// 
// This free software is complimented by the author to accompany the textbook.
// E-mail: altacz@gmail.com.

//                                                                                   
// DESCRIPTION: A Subroutine to find the roots of a quadratic equation of            
//   the form : x*x + p * x + q = 0.                                                 
//                                                                                   
// ON ENTRY                                                                          
//   p, q :: Coefficients of the quadratic equation;                                 
//                                                                                   
// ON EXIT                                                                           
//   re   :: Array of length 2 containing real parts of the roots: re1, re2;         
//   im   :: Array of length 2 containing imaginary parts of the roots: im1, im2.    
//                                                                                   
// USES                                                                              
//   SQRT :: Built-in Intrinsic function to evaluate the square root of a real value.
//                                                                                   
// REVISION DATE :: 03/18/2024                                                       
//                                                                                   
// ==================================================================================
    double d = p * p - 4.0 * q;
    if (d < 0.0) {
        d = std::sqrt(-d);
        re1 = -p / 2.0;
        re2 = re1;
        im1 = -d / 2.0;
        im2 = -im1;
    } else {
        d = std::sqrt(d);
        re1 = (-p - d) / 2.0;
        im1 = 0.0;
        re2 = (-p + d) / 2.0;
        im2 = 0.0;
    }
}

// =============================================================================
//  The main program to test "quadratic_eq"
//  ============================================================================
int main() {
    double p, q, xr1, xi1, xr2, xi2;
    std::cout << "Type in values for p and q" << std::endl;
    std::cin >> p >> q;
    quadratic_eq(p, q, xr1, xi1, xr2, xi2);
    std::cout << "1st Root " << xr1 << " + " << xi1 << " i" << std::endl;
    std::cout << "2nd Root " << xr2 << " + " << xi2 << " i" << std::endl;
    return 0;
}

 