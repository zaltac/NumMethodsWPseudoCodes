#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

double exact(double x) {
// ==============================================================================
// DESCRIPTION: A function subprogram providing the true solution y=f(x) for testing the module. 
//
// ARGUMENTS:
//      x   :: A real input, independent variable.
// ==============================================================================
    return x * x * (x * x - 0.5);
}

void coeffs(double x, double &p, double &q, double &r, double &f) {
// ==============================================================================
// DESCRIPTION: A user-defined module suppling the coefficients p(x), q(x), r(x) and
//    rhs g(x) of the linear ODE given in the following form: 
//         p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]
//
// ON ENTRY
//    x   :: Independent variable (a <=x<= b);
//
// ON EXIT
//   p, q, r, abd g :: Coefficients & rhs evaluated at x. 
//
// REVISION DATE :: 03/08/2025
// ==============================================================================
    p = x * x;
    q = -5.0 * x;
    r = 8.0;
    f = 0.0;
}

void tridiagonal(int s1, int sn, std::vector<double> &b, std::vector<double> &d, std::vector<double> &a, std::vector<double> &c, std::vector<double> &x) {
// ==================================================================================
// CODE2.13-tridiagonal.cpp. A C++ module implementing Pseudocode 2.13.                           
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
// DESCRIPTION: A module to solve a tridiagonal system of linear equations using Thomas algorithm.                           
//                                                                                             
// ON ENTRY                                                                                    
//    s1 :: Subscript of the first unknown (usually 1);                                        
//    sn :: Subscript of the last unknown (usually no. of eqs, n);                             
//     b :: Array of length n containing coefficients of below diagonal elements;              
//     d :: Array of length n containing coefficients of diagonal elements;                    
//     a :: Array of length n containing coefficients of above diagonal elements;              
//     c :: Array of length n containing coefficients of rhs.                                  
//                                                                                             
// ON EXIT                                                                                     
//     x :: An array of length n containing the solution.                                      
//                                                                                             
// REVISION DATE :: 03/18/2024                                                                 
// ==================================================================================
    for (int i = s1 + 1; i <= sn; i++) {
        double ratio = b[i] / d[i - 1];
        d[i] = d[i] - ratio * a[i - 1];
        c[i] = c[i] - ratio * c[i - 1];
    }

    x[sn] = c[sn] / d[sn];
    for (int i = sn - 1; i >= s1; i--) {
        x[i] = (c[i] - a[i] * x[i + 1]) / d[i];
    }
}

std::vector<double> LBVP_SOLVE(int neq, std::vector<double> &x, std::vector<double> &alpha, std::vector<double> &beta, std::vector<double> &gamma) {
// ==================================================================================
// CODE10.1-LBVP_SOLVE.cpp. A C++ module implementing Pseudocode 10.1.                            
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
// DESCRIPTION: A module to find approximate solution of the following linear                  
//   differential equation using the Finite Difference Method:                                 
//         p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]                                   
//   subject to                                                                                
//         alpha1 * y'(a)+ beta1 * y(a) = gamma1                                               
//         alpha2 * y'(b)+ beta2 * y(b) = gamma2                                               
//                                                                                             
// CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1.  
//                                                                                             
// ON ENTRY                                                                                    
//   neq  :: Number of (equations) grid poinds;                                                
//    x   :: Array of length neq containing the abscissa of the grid points;                   
//   alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                   
//          the boundary conditions as stated above;                                           
//                                                                                             
// ON EXIT                                                                                     
//    y   :: Array of length neq containing the approximate solution.                          
//                                                                                             
// USES                                                                                        
//  COEFFS  :: A module containing the coefficients and rhs of the linear ODE, i.e., p(x), q(x), r(x), g(x).
//  TRIDIAGONAL:: A module solving a tridiagonal system of equations using the Thomas algorithm
//                                                                                             
// REVISION DATE :: 03/08/2025                                                                 
// ==================================================================================    
    int nbc1 = static_cast<int>(alpha[0]);
    int nbc2 = static_cast<int>(alpha[1]);
    std::vector<double> bb(neq);
    std::vector<double> aa(neq);
    std::vector<double> dd(neq);
    std::vector<double> cy(neq);
    double h = x[1] - x[0];
    double h2 = h * h;

    for (int k = 0; k < neq; k++) {
        double pxk, qxk, rxk, fxk;
        coeffs(x[k], pxk, qxk, rxk, fxk);
        dd[k] = rxk - 2.0 * pxk / h2;
        aa[k] = pxk / h2 + 0.5 * qxk / h;
        bb[k] = pxk / h2 - 0.5 * qxk / h;
        cy[k] = fxk;
    }

    if (nbc1 == 0) {
        dd[0] = 1.0;
        aa[0] = 0.0;
        bb[0] = 0.0;
        cy[0] = gamma[0] / beta[0];
    } else {
        aa[0] = aa[0] + bb[0];
        dd[0] = dd[0] + 2.0 * h * bb[0] * beta[0] / alpha[0];
        cy[0] = cy[0] + 2.0 * h * bb[0] * gamma[0] / alpha[0];
    }

    if (nbc2 == 0) {
        dd[neq - 1] = 1.0;
        aa[neq - 1] = 0.0;
        bb[neq - 1] = 0.0;
        cy[neq - 1] = gamma[1] / beta[1];
    } else {
        bb[neq - 1] = aa[neq - 1] + bb[neq - 1];
        dd[neq - 1] = dd[neq - 1] - 2.0 * h * aa[neq - 1] * beta[1] / alpha[1];
        cy[neq - 1] = cy[neq - 1] - 2.0 * h * aa[neq - 1] * gamma[1] / alpha[1];
    }

    std::vector<double> y(neq);
    tridiagonal(0, neq - 1, bb, dd, aa, cy, y);

    return y;
}

int main() {
// ==============================================================================
//  The main program to test LBVP_SOLVE.CPP
// ==============================================================================
    int n;
    std::cout << "Enter no of intervals: ";
    std::cin >> n;
    int neq = n + 1; // no. of equations

    std::vector<double> x(neq);
    std::vector<double> y;
    double xa = 1.0, xb = 2.0;
    double h = (xb - xa) / n;
    for (int i = 0; i < neq; i++) {
        x[i] = xa + i * h;
    }

    std::vector<double> alpha = {1.0, 1.0};
    std::vector<double> beta  = {2.0, -2.0};
    std::vector<double> gamma = {4.0, 2.0};

    y = LBVP_SOLVE(neq, x, alpha, beta, gamma);

    std::cout << std::setw(15) << "      x" << std::setw(10) << "Exact" << std::setw(14) << "N.Approx" << std::setw(16) << "Abs Error" << std::endl;
    for (int i = 0; i < neq; i++) {
        double error = std::fabs(y[i] - exact(x[i]));
        std::cout << std::setw(10) << std::fixed << std::setprecision(4) << x[i]
                  << std::setw(12) << std::fixed << std::setprecision(7) << exact(x[i])
                  << std::setw(12) << std::fixed << std::setprecision(7) << y[i]
                  << std::setw(16) << std::scientific << std::setprecision(5) << error << std::endl;
    }

    return 0;
}