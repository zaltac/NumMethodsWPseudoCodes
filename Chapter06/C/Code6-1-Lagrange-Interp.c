#include <stdio.h>

double LAGRANGE_P(int n, double xval, double x[], double L[]) {
// ==============================================================================
// CODE6.1-LAGRANGE_P.C. A C-module implementing Pseudocode 6.1.                                
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
// DESCRIPTION: A C module to compute the Lagrange polynomials at x=xval.                        
//                                                                                             
// ON ENTRY                                                                                    
//   n    :: The number of data in the set minus 1;                                            
//   xval :: x value at which the dataset is to be interpolated;                               
//   x    :: Array of length n+1 containing abscissas, k=0,1,2,...,n.                          
//                                                                                             
// ON EXIT                                                                                     
//   L    :: Array of length (n+1) containing Lagrange polynomials                             
//           that is, L(k) = L_k(xval) for k=0,1,2,...,n.                                      
//                                                                                             
// REVISION DATE :: 02/28/2025  
// ==============================================================================
    for (int j = 0; j <= n; j++) {
        L[j] = 1.0;
    }

    for (int k = 0; k <= n; k++) {
        for (int i = 0; i <= n; i++) {
            if (i != k) {
                L[k] *= (xval - x[i]) / (x[k] - x[i]);
            }
        }
    }

    return L[n+1];
}

double LAGRANGE_EVAL(int n, double xval, double x[], double f[]) {
// ==============================================================================
// CODE6-1-LAGRANGE_EVAL.C. A C-module for implementing Pseudocode 6.1.
//     
// NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
// First Edition. (c) By Zekeriya ALTAÇ (2024).
// ISBN: 9781032754741 (hbk)
// ISBN: 9781032756424 (pbk)
// ISBN: 9781003474944 (ebk)
//
// DOI : 10.1201/9781003474944
// C&H/CRC PRESS, Boca Raton & London. 
//  
// This free software is complimented by the author to accompany the textbook.
// E-mail: altacz@gmail.com
//
// DESCRIPTION: A function to evaluate Lagrange interpolating polynomial for an arbitrary xval, f=f(xval)=fval 
//  within x0 <= xval <= xn.
//
// ARGUMENTS
//   n    :: The number of data in the set minus 1;
//   xval :: x value at which the dataset is to be interpolated;
//   x    :: Array of length (n+1) containing abscissas, k=0,1,2,...,n;
//   f    :: Array of length (n+1) containing ordinates, k=0,1,2,...,n.
//
// USES
//   LAGRANGE_P :: A C-module for generating the Lagrange polynomials.
//
// REVISION DATE :: 02/28/2025
// ==============================================================================
    double L[n + 1];

    LAGRANGE_P(n, xval, x, L); 
    
    double fval = 0.0;

    for (int k = 0; k <= n; k++) {
        fval += L[k] * f[k];
    }

    return fval;
}


// ==============================================================================
//  The main program to test subprograms LAGRANGE_EVAL.C 
// ==============================================================================

int main() {
    int n = 3;
    double x[] = {0.08, 0.25, 0.5, 0.9};
    double f[] = {0.25, 0.625, 0.81, 0.43};
    double xval;

    printf("Enter x: ");
    scanf("%lf", &xval);

    printf("fval = %lf\n", LAGRANGE_EVAL(n, xval, x, f));

    return 0;
}