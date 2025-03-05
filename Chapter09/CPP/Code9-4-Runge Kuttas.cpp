#include <iostream>
#include <cmath>
#include <iomanip>


double FCN(double x, double y) {
// ==============================================================================
// DESCRIPTION: A function subprogram providing y'=f(x,y)
//
// ARGUMENTS:
//      x, y  :: Double input values.
//
// ==============================================================================
    return -y / ( x + 1.0 );
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
    return 2.0 /( x + 1.0 );
}


void DRV_RK(int n, double h, double x0, double y0, double &x, double &y) {
// ==================================================================================
// CODE9.4-DRV_RK.CPP. A C++ module implementing driver module of Pseudocode 9.4.                                  
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
// DESCRIPTION: A driver module employing one-step RK2, RK3, or RK4 scheme.               
//                                                                                             
// ON ENTRY                                                                                    
//  n     :: Order of Runge-Kutta scheme;                                                      
//  h     :: Step size (it must be uniform);                                                   
//  x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
//                                                                                             
// ON EXIT                                                                                     
//  x,y   :: Current estimates, x^(p+1) and y^(p+1).                                           
//                                                                                             
// USES                                                                                                
//   FCN  :: User-defined external function providing y'=f(x,y).                               
//                                                                                             
// REVISION DATE :: 03/05/2025                                                                 
// ==================================================================================
    double x1, xh, ym, xk, xk1, xk2, xk3, xk4;
    double hlf = 0.50;

    xh = x0 + 0.5 * h;
    x1 = x0 + h;

    switch (n) {
        case 2:  // Case of RK2
            xk1 = h * FCN(x0, y0);
            ym = y0 + xk1;
            xk2 = h * FCN(x1, ym);
            xk = hlf * (xk1 + xk2);
            break;

        case 3:  // Case of RK3
            xk1 = h * FCN(x0, y0);
            ym = y0 + hlf * xk1;
            xk2 = h * FCN(xh, ym);
            ym = y0 - xk1 + 2.0 * xk2;
            xk3 = h * FCN(x1, ym);
            xk = (xk1 + 4.0 * xk2 + xk3) / 6.0;
            break;

        case 4:  // Case of RK4
            xk1 = h * FCN(x0, y0);
            ym = y0 + hlf * xk1;
            xk2 = h * FCN(xh, ym);
            ym = y0 + hlf * xk2;
            xk3 = h * FCN(xh, ym);
            ym = y0 + xk3;
            xk4 = h * FCN(x1, ym);
            xk = (xk1 + 2.0 * xk2 + 2.0 * xk3 + xk4) / 6.0;
            break;

        default:
            std::cerr << "PROGRAM DOES NOT HANDLE CASE OF N=" << n << std::endl;
            exit(1);
    }

    y = y0 + xk;
    x = x1;
}



void Runge_Kutta(int n, double h, double x0, double y0, double xlast) {
// ==================================================================================
// CODE9.4-Runge_Kutta.cpp. A C++ module implementing Pseudocode 9.4.                             
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
// DESCRIPTION: A module to estimate the solution of a first order IVP on [x0,xlast]       
//   using the 2nd to 4th order Runge-Kutta scheme. Numerical estimates are printed out, not stored.           
//                                                                                             
// ON ENTRY                                                                                    
//    n    :: Order of the Runge-Kutta scheme;                                                 
//    h    :: Step size (it must be uniform);                                                  
//   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps
//  xlast  :: End point of the solution interval.                                              
//                                                                                             
// Other Internal Variables                                                                    
//   x, y  :: Current estimates, x^(p+1) and y^(p+1).                                          
//                                                                                             
//  USES                                                                                       
//    abs  :: Built-in Intrinsic function returning the absolute value of a real value.        
//   DRV_RK:: A driver subprogram performing one-step RK scheme.                               
//                                                                                             
// REVISION DATE :: 03/05/2025                                                                 
// ==================================================================================
    double yt, aerr, x, y;

    std::cout << std::setw(4) << "x0" << std::setw(12) << "y0" << std::endl;
    std::cout << std::setw(4) << x0 << std::setw(12) << y0 << std::endl;

    x = x0;
    while (x < xlast) {
        DRV_RK(n, h, x0, y0, x, y);
        yt = exact(x);
        aerr = std::abs(y - yt);
        std::cout << std::setw(12) << x << std::setw(12) << yt << std::setw(12) << y << std::setw(12) << aerr << std::endl;
        x0 = x;
        y0 = y;
    }
}


int main() {
// ==============================================================================
//  The main program to test Runge_Kutta.CPP
// ==============================================================================
    int n;
    double x0 = 0.0;
    double xlast = 1;
    double y0 = 2.0;
    double h = 0.1;

    std::cout << "Enter Order of RK method: ";
    std::cin >> n;

    Runge_Kutta(n, h, x0, y0, xlast);

    return 0;
}
