#include <math.h>
#include <stdio.h>

//  ==================================================================================
//  USER-DEFINED FUNCTION "FUNC" OF ONE-VARIABLE
//  ==================================================================================
double FUNC(double x) {
    return 25000.0 / (-57.0 + x) - 5.2e6 / (x * x);
}


// ==================================================================================
// CODE5.2-RICHARDSON.C. A C module implementing Pseudocode 5.2.                                  
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
// DESCRIPTION: A C module to compute the first derivative of an explicitly                      
//     defined function using Richardson's extrapolation.                                    
//                                                                                             
// ON ENTRY                                                                                    
//     x0  :: Point at which derivative is to be computed;                                     
//     h   :: Initial interval size;                                                           
//     eps :: Tolerance desired.                                                               
//                                                                                             
// ON EXIT                                                                                     
//     D   :: A matrix containing the Richardson's table (0..n, 0..n)                        
//     nr  :: Size of the table;                                                               
//   deriv :: Estimated derivative.                                                            
//                                                                                             
// USES                                                                                        
//    fabs :: Built-in Intrinsic function returning the absolute value of a real value.        
//                                                                                             
// ALSO REQUIRED                                                                               
//    FUNC  :: User-defined external function providing the nonlinear equation.                
//                                                                                             
// REVISION DATE :: 06/13/2024                                                                 
// ==================================================================================
void Richardson(double x0, double h, double eps, double D[11][11], int* nr, double* deriv) {
    int k = 0, m;
    double err = 1.0;

    while (err > eps) {
        D[k][0] = (FUNC(x0 + h) - FUNC(x0 - h)) / (2.0 * h);  // 1st derivative
        for (m = 1; m <= k; m++) {
            D[k][m] = (pow(4, m) * D[k][m - 1] - D[k - 1][m - 1]) / (pow(4, m) - 1);
        }
        if (k >= 1) {  // Estimate diagonalwise differentiation error
            err = fabs(D[k][k] - D[k - 1][k - 1]);
            printf("D(%d,%d)-D(%d,%d)=%lf\n", k, k, k - 1, k - 1, err);
        }
        h /= 2;
        k++;
    }

    *nr = k - 1;
    *deriv = D[*nr][*nr];
}


// ==============================================================================
//  The main program to test module Richardson
// ==============================================================================
int main() {
    double D[11][11], x0 = 150.0, h = 5.0, err = 1.0, eps = 1.0e-6, deriv;
    int nr, k, m;

    Richardson(x0, h, eps, D, &nr, &deriv);

    for (k = 0; k <= nr; k++) {
        printf("%d ", k);
        for (m = 0; m <= k; m++) {
            printf("%lf ", D[k][m]);
        }
        printf("\n");
    }

    printf("---------------------------------\n");
    printf("Derivative is=%lf\n", deriv);
    printf("---------------------------------\n");

    return 0;
}