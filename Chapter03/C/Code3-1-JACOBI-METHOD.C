#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// ==================================================================================
// CODE3.1-ENORM.C. A C module implementing ENORM of Pseudocode 3.1.                                      
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
// DESCRIPTION: A function module to compute Euclidean (L2) norm of a vector.                  
//                                                                                             
// ARGUMENTS                                                                                   
//     n  :: The length of an input vector;                                                    
//     x  :: A vector (array) of length n.                                                     
//                                                                                             
// USES                                                                                        
//   SQRT :: Built-in Intrinsic function returning the square root of a real value.            
//                                                                                             
// REVISION DATE :: 11/09/2024                                                                  
// ==================================================================================
double ENorm(int n, double *x) {
    double sums;
    sums=0.0;
    for (int i = 0; i <= n; i++) {
        sums += x[i] * x[i];
    }
    sums = sqrt(sums); 
    return sums;
}


// ==================================================================================
// CODE3.1-JACOBI.C. A C module implementing JACOBI_DRV module in Pseudocode 3.1.                                      
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
// DESCRIPTION: A module to to perform one step Jacobi iteration and compute               
//    the Euclidean norm of the displacement vector.                                           
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Number of equations (size of A);                                                  
//    A   :: Input coefficient matrix (nxn);                                                  
//    b   :: Array of length n containing the right-hand;                                      
//    x   :: Array of length n containing the estimate at (p+1)'th step;                       
//    xo  :: Array of length n containing the estimate at p'th step.                           
//                                                                                             
// ON EXIT                                                                                     
//    x   :: Array of length n containing the estimated solution;                              
//    del :: Maximum absolute error achieved.                                                  
//                                                                                             
// USES                                                                                        
//  ENORM:: User-defined function calculating the Euclidean vector (L2 norm) of a vector.      
//                                                                                             
// REVISION DATE :: 11/09/2024                                                                  
// ==================================================================================
void jacobi_drv(int n, double **a, double *b, double *xo, double *x, double *delta) {
    double *dx = (double *)calloc(n, sizeof(double));  
    for (int i = 0; i < n; i++) {
        double sums = 0.0;
        for (int k = 0; k < n; k++) {
            if (i != k) {
                sums += a[i][k] * xo[k];
            }
        }
        x[i] = (b[i] - sums) / a[i][i];
        dx[i] = x[i] - xo[i];
    }
    for (int k = 0; k < n; k++) {
        dx[k] = x[k] - xo[k];
    }
    *delta = ENorm(n,dx);
    free(dx);
}


// ==================================================================================
// CODE3.1-JACOBI.C. A C module implementing Pseudocode 3.1.                                      
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
// DESCRIPTION: A module to iteratively solve Ax=b using the Jacobi method.                
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Number of equations (size of A);                                                  
//    A   :: Input coefficient matrix (nxn);                                                  
//    b   :: Array of length n containing the right-hand;                                      
//    x   :: Array of length n containing the estimate at (p+1)'th step;                       
//    xo  :: Array of length n containing the initial guess, or iterates at estimate at p'th st
//   eps  :: Convergence tolerance;                                                            
//  maxit :: Maximum permitted number of iterations.                                           
//                                                                                             
// ON EXIT                                                                                     
//    x   :: Array of length n containing the estimated solution;                              
//  iter  :: Total number of iterations performed;                                             
//  error :: L2 norm of the displacement vector.                                               
//                                                                                             
// USES                                                                                        
//   JACOBI_DRV :: Accompanying module performing one step Jacobi iteration.               
//                                                                                             
// REVISION DATE :: 11/09/2024                                                                  
// ==================================================================================
void jacobi(int n, double eps, double **a, double *b, double *xo, int maxit, int iprint) {
    double del1, del0 = 1.0;
    int p = 0;
    double *x = (double *)calloc(n, sizeof(double));
    while (del0 > eps && p<maxit) {
        p++;
        jacobi_drv(n, a, b, xo, x,&del1);
        if (iprint == 1) {
            printf("iter=%d  Error= %10.7f  Ratio= %f\n", p, del1, del1 / del0);
        }
        for(int i=0; i<n; i++) {
            xo[i] = x[i];
        }
        del0 = del1;
    }
    if (p == maxit) {
        printf("\nJacobi method failed to converge after %d iterations\nwithin the specified EPS tolerance.\n", maxit);
    }
    printf("\n------- Solution ---------------\n");
    for (int i = 0; i < n; i++) {
        printf(" x(%d) = %f\n", i + 1, x[i]);
    }
    printf("\n-----------------------\n");
    printf("Total no of iterations = %d Maximum Error = %f\n", p, del0);
    free(x);
}

// ==============================================================================
//  The main program to test JACOBI.C
// ==============================================================================
int main() {
    int n = 10;
    double **a = (double **)calloc(n, sizeof(double*));
    for(int i=0; i<n; i++) {
        a[i] = (double *)calloc(n, sizeof(double));
    }
    double *b = (double *)calloc(n, sizeof(double));
    double *xo = (double *)calloc(n, sizeof(double));
    int iprint = 1;

    a[0][0] = 2.0; a[0][1] = -1.0;
    a[1][0] = -1.0; a[1][1] = 2.0; a[1][2] = -1.0;
    a[2][1] = -1.0; a[2][2] = 2.0; a[2][3] = -1.0;
    a[3][2] = -1.0; a[3][3] = 2.0; a[3][4] = -1.0;
    a[4][3] = -1.0; a[4][4] = 2.0; a[4][5] = -1.0;
    a[5][4] = -1.0; a[5][5] = 2.0; a[5][6] = -1.0;
    a[6][5] = -1.0; a[6][6] = 2.0; a[6][7] = -1.0;
    a[7][6] = -1.0; a[7][7] = 2.0; a[7][8] = -1.0;
    a[8][7] = -1.0; a[8][8] = 2.0; a[8][9] = -1.0;
    a[9][8] = -1.0; a[9][9] = 2.0;
    b[9] = 2.2;

    printf("* Coefficient Matrix *" " * RHS *" " * Initial Guess *\n");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%6.1f ", a[i][j]);
        }
        printf("%12.1f %12.1f\n", b[i], xo[i]);
    }    

    double eps = 1.0e-5;
    int maxit = 999;
 
    jacobi(n, eps, a, b, xo, maxit, iprint);

    for(int i=0; i<n; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(xo);

    return 0;
}