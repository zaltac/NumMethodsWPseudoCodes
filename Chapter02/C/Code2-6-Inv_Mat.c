#include <stdio.h>
#include <stdlib.h>
 
// ==================================================================================
// CODE2.6-INV_MAT.C. A C module implementing Pseudocode 2.6.                                     
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
// DESCRIPTION: A module to find the inverse of a square matrix (with no pivoting).        
//                                                                                             
// ON ENTRY                                                                                    
//    n  :: Dimension attribute of input matrix A;                                             
//    A  :: An input matrix (nxn).                                                             
//                                                                                             
// ON RETURN                                                                                    
//    B  :: Inverse of A (nxn).                                                                
//                                                                                             
// REVISION DATE :: 03/18/2024                                                                 
// ==================================================================================
double** inv_mat(double** a, int n) {
    double **b = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        b[i] = (double*)malloc(n * sizeof(double));
    }

    // Initialize the inverse matrix as the identity matrix
    for (int i = 0; i < n; i++) {
        b[i][i] = 1.0;
    }

    for (int j = 0; j < n; j++) {
        double p = 1.0 / a[j][j];
        for (int k = 0; k < n; k++) {
            b[j][k] *= p;
            a[j][k] *= p;
        }

        for (int i = 0; i < n; i++) {
            if (i != j) {
                double s = a[i][j];
                for (int k = 0; k < n; k++) {
                    b[i][k] -= s * b[j][k];
                    a[i][k] -= s * a[j][k];
                }
            }
        }
    }

    return b;
}

// ==============================================================================
//  The main program to test Inv_Mat.c
// ==============================================================================
int main() {
    int n = 5;
    double **a = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        a[i] = (double*)malloc(n * sizeof(double));
    }
    a[0][0] = 1;  a[0][1] = 2; a[0][2] = 1; a[0][3] = 2; a[0][4] = 3;
    a[1][0] = 11; a[1][1] = -1; a[1][2] = 1; a[1][3] = 4; a[1][4] = 1;
    a[2][0] = 4;  a[2][1] = -1; a[2][2] = 1; a[2][3] = 1; a[2][4] = -1;
    a[3][0] = -3; a[3][1] = 1; a[3][2] = -8; a[3][3] = -1; a[3][4] = 5;
    a[4][0] = -1; a[4][1] = 1; a[4][2] = 1; a[4][3] = 1; a[4][4] = 1;

    printf("\n ********* Input Matrix *********\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.5f ", a[i][j]);
        }
        printf("\n");
    }

    double** ainv = inv_mat(a, n);

    printf("\n -------- matrix A(-1) --------\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.5f ", ainv[i][j]);
        }
        printf("\n");
    }
    printf("-------------------------------\n");


    return 0;
}