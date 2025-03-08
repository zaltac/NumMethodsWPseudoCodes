#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double exact(double x) {
// ==============================================================================
// DESCRIPTION: A function module providing the true solution y=f(x) for testing the module. 
//
// ARGUMENTS:
//      x   :: A real input, independent variable.
// ==============================================================================
    return x * x * (x * x - 0.5);
}

void coeffs(double x, double *p, double *q, double *r, double *f) {
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
    *p = x * x;
    *q = -5.0 * x;
    *r = 8.0;
    *f = 0.0;
}

void tridiagonal(int s1, int sn, double *b, double *d, double *a, double *c, double *x) {
// ==================================================================================
// CODE2.13-tridiagonal.c. A C module implementing Pseudocode 2.13.                 
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

double* LBVP_SOLVE(int neq, double *x, double *alpha, double *beta, double *gamma) {
// ==================================================================================
// CODE10.1-LBVP_SOLVE.c. A C module implementing Pseudocode 10.1.                                
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
// DESCRIPTION: A function module to find approximate solution of the following linear                  
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
    int nbc1 = (int)alpha[0];
    int nbc2 = (int)alpha[1];
    double *bb = (double *)calloc(neq, sizeof(double));
    double *aa = (double *)calloc(neq, sizeof(double));
    double *dd = (double *)calloc(neq, sizeof(double));
    double *cy = (double *)calloc(neq, sizeof(double));
    double h = x[1] - x[0];
    double h2 = h * h;

    for (int k = 0; k < neq; k++) {
        double pxk, qxk, rxk, fxk;
        coeffs(x[k], &pxk, &qxk, &rxk, &fxk);
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

    double *y = (double *)calloc(neq, sizeof(double));
    tridiagonal(0, neq - 1, bb, dd, aa, cy, y);

    free(bb);
    free(aa);
    free(dd);
    free(cy);
    
    return y;
}



int main() {
// ==============================================================================
//  The main program to test LBVP_SOLVE.C
// ==============================================================================
    int n;
    printf("Enter no of intervals: ");
    scanf("%d", &n);
    int neq = n + 1; // no. of equations

    double *x = (double *)calloc(neq, sizeof(double));
    double *y;
    double xa = 1.0, xb = 2.0;
    double h = (xb - xa) / n;
    for (int i = 0; i < neq; i++) {
        x[i] = xa + i * h;
    }

    double alpha[2] = {1.0, 1.0};
    double beta[2] = {2.0, -2.0};
    double gamma[2] = {4.0, 2.0};

    y = LBVP_SOLVE(neq, x, alpha, beta, gamma);

    printf("%15s%10s%14s%16s\n", "      x", "Exact", "N.Approx", "Abs Error");
    for (int i = 0; i < neq; i++) {
        double error = fabs(y[i] - exact(x[i]));
        printf("%10.4f%12.7f%12.7f%16.5e\n", x[i], exact(x[i]), y[i], error);
    }

    free(x);
    free(y);
    return 0;
}