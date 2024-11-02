#include <stdio.h>
#include <stdlib.h>

double* Ax(int n, double** A, double* x) {
// ==================================================================================
// CODE2.5-Ax.C. A C module implementing Pseudocode 2.5.                                          
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
// DESCRIPTION: A module to perform A * x = b matrix-vector multiplication.                
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Dimension attributes of input/output matrices;                                    
//    A   :: An input matrix of size nxn;                                                     
//    x   :: An input vector of length n.                                                      
//                                                                                             
// ON RETURN                                                                                     
//    b   :: The output vector of length n.                                                    
//                                                                                             
// REVISION DATE :: 03/18/2024                                                                 
// ==================================================================================
    double* b = (double*)calloc(n, sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * x[j];
        }
    }
    return b;
}


// ==============================================================================
//  The main program to test Ax.c
// ==============================================================================
int main() {
    int n = 4;
    double** A = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        A[i] = (double*)malloc(n * sizeof(double));
    }
    A[0][0] = 1.0; A[0][1] = 4.0; A[0][2] = 3.0; A[0][3] = -2.0;
    A[1][0] = 2.0; A[1][1] = 1.0; A[1][2] = 2.0; A[1][3] = 3.0;
    A[2][0] = 3.0; A[2][1] = 2.0; A[2][2] = 1.0; A[2][3] = 4.0;
    A[3][0] = -2.0; A[3][1] = 3.0; A[3][2] = 2.0; A[3][3] = 1.0;
    double x[] = {2.0, 5.0, 3.0, -3.0};

    printf(" Input Matrix A\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n Input Vector x\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", x[i]);
    }
    printf("\n");

    double* b = Ax(n, A, x);

    printf("\n ------ A(n,n)X(n) product is \n");
    printf("\n Output Vector b\n");

    for (int i = 0; i < n; i++) {
        printf("%f ", b[i]);
    }
    printf("\n");

    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);

    return 0;
}
