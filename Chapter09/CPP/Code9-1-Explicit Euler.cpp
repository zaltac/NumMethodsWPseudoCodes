#include <iostream>
#include <cmath>
 


double fcn(double x, double y) {
// ==============================================================================
// DESCRIPTION: A function subprogram providing y'=f(x,y)
//
// ARGUMENTS:
//      x, y  :: Real input values.
//
// ==============================================================================
    return -y/(x + 1.0);
}

double exact(double x) {
// ==============================================================================
// DESCRIPTION: A function subprogram providing the true solution y=(x) for 
//    testing the module. 
//
// ARGUMENTS:
//      x   :: A double input, independent variable.
//
// ==============================================================================
    return 2.0 / (x + 1.0 );
}


void Explicit_Euler(double h, double x0, double y0, double xlast) {
// ==================================================================================
// CODE9.1-Explicit_Euler.cpp. A C++ module implementing Pseudocode 9.1.                          
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
// DESCRIPTION: A C++ module to estimate the solution of a first order IVP on [x0,xlast]           
//   using the Explicit Euler method. Numerical estimates are printed out, not stored.         
//   With minor modifications, this PROGRAM can also be used to solve explicit methods         
//   such as MIDPOINT RULE and MODIFIED EULER.                                                 
//                                                                                             
// ON ENTRY                                                                                    
//  h     :: Step size (it must be uniform);                                                   
//  x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
//  xlast :: End point of the solution interval.                                               
//                                                                                             
// Other Variables                                                                             
//  x,y   :: Current estimates, x^(p+1) and y^(p+1).                                           
//                                                                                             
// USES                                                                                        
//   abs  :: Built-in intrinsic function returning the absolute value of a real value;
//   FCN  :: User-defined external function providing y'=f(x,y).                               
//                                                                                             
// REVISION DATE :: 03/05/2024                                                                 
// ==================================================================================
    double x = x0;
    double y = y0;

    std::cout << "    x          y          Error" << std::endl;
    std::cout << "---------------------------------" << std::endl;
    std::cout << x0 << "    " << y0 << std::endl;

    while (x <= xlast) {
        x = x0 + h;
        y = y0 + h * fcn(x0, y0);  // Explicit-Euler
        std::cout << x << "    " << y << "    " << std::abs(y - exact(x)) << std::endl;
        x0 = x;
        y0 = y;
    }
}



int main() {
// ==============================================================================
//  The main program to test Explicit_Euler.CPP
// ==============================================================================
    double x0 = 0.0;
    double y0 = 2.0;
    double xlast = 1.0;
    double h = 0.10;

    Explicit_Euler(h, x0, y0, xlast);

    return 0;
}