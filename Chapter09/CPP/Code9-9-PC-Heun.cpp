#include <iostream>
#include <cmath>
#include <iomanip>

double FCN(double x, double y) {
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



void PC_Heun(double h, double x0, double y0, double xlast) {
// ==================================================================================
// CODE9.9-PC_HEUN.CPP. A C++ module implementing Pseudocode 9.9.                         
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
// DESCRIPTION: A CPP module to estimate the solution of a first order IVP on [x0,xlast]           
//   using the Heun's Predictor-Correcter method. Numerical estimates are printed out, not stored.                                              
//                                                                                             
//  ON ENTRY                                                                                   
//   h     :: Step size (it must be uniform);                                                  
//   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps
//   xlast :: End point of the solution interval.                                              
//                                                                                             
//  Other Internal Variables                                                                   
//   x,y   :: Current estimates, x^(p+1) and y^(p+1).                                          
//                                                                                             
// USES                                                                                        
//   abs  :: Built-in CMATH library function returning the absolute value of a real value.         
//   FCN  :: User-defined external function providing y'=f(x,y).                               
//                                                                                             
// REVISION DATE :: 03/07/2025                                                                 
// ==================================================================================
    double x = x0;
    double y = y0;
    double ys, k1, k2, err;

    std::cout << std::fixed << std::setprecision(7) << std::setw(12) << x0 << "  " << std::setw(12) << y0 << "  " << std::scientific << std::setw(14) << std::abs(y0 - exact(x0)) << std::endl;

    while (x < xlast) {
        // Predictor step
        k1 = h * FCN(x0, y0);
        ys = y0 + k1;

        // Corrector step
        x = x0 + h;
        k2 = h * FCN(x, ys);
        y = y0 + 0.50 * (k1 + k2);
        err = std::abs(y - exact(x));

        std::cout << std::fixed << std::setprecision(7) << std::setw(12) << x << "  " << std::setw(12) << y << "  " << std::scientific << std::setw(14) << err << std::endl;

        x0 = x;
        y0 = y;
    }
}



int main() {
// ==============================================================================
//  The main program to test PC_Heun.CPP
// ==============================================================================
    double x0 = 0.0;
    double y0 = 2.0;
    double h = 0.1;
    double xlast = 1.0;

    PC_Heun(h, x0, y0, xlast);

    return 0;
}

