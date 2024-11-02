#include <stdio.h>
#include <stdlib.h>

// ==================================================================================
// CODE2.4-MAT_MUL.C. A C module implementing Pseudocode 2.4.                                     
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
//  DESCRIPTION: A module to find A*B=C matrix multiplication.                             
//                                                                                             
//  ON ENTRY                                                                                   
//   m,p,n :: Dimension attributes of input/output matrices;                                   
//      A  :: An input matrix of size mxp;                                                    
//      B  :: An input matrix of size pxn.                                                    
//                                                                                             
//  ON RETURN                                                                                    
//      C  :: The output matrix of size mxn.                                                  
//                                                                                             
//  REVISION DATE :: 03/18/2024                                                                
// ==================================================================================
double** mat_mul(double** A, double** B, int m, int pA, int pB, int n) {
    double** C = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        C[i] = (double*)malloc(n * sizeof(double));
    }
    printf("\nInside mat_mul \n");
    printf("\n%d",m);
    if (pA == pB) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < pA; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
                printf("C(%d)=%f\n", i, C[i][j]);                
            }
        }
    }
    else {
        printf("Sorry, cannot multiply A and B.\n");
        printf("Execution halted!!!\n");
        exit(0);
    }

    return C;
}

// ==============================================================================
//  The main program to test mat_mul.c
// ==============================================================================
int main() {
    double A[2][4] = {{1.0, -3.0, 5.0, 1.0}, {-2.0, 4.0, 1.0, 2.0}};
    double B[4][3] = {{3.0, 1.0, -2.0}, {2.0, 0.0, 4.0}, {-1.0, 1.0, -3.0}, {2.0, 5.0, 3.0}};

    printf("Input Matrix A \n\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }

    printf("\nInput Matrix B \n\n");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", B[i][j]);
        }
        printf("\n");
    }
    printf("\nEnter mat_mul \n");
    
    double** C = mat_mul((double**)A, (double**)B, 2, 4, 4, 3);

    printf("\n------ A(m,p)xB(p,n) matrix product is ------\n\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", C[i][j]);
        }
        printf("\n");
    }


    return 0;
}