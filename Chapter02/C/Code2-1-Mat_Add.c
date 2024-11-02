#include <stdio.h>

void mat_add(int m, int n, double A[m][n], double B[m][n], double C[m][n]);

// ==============================================================================
//  The main program to test MAT_ADD.C
// ==============================================================================
const int n = 3;
const int m = 3;
    
int main() {
    double C[m][n];
    int i, j;

    // Initialize matrices A and B
    double A[3][3] ={ {1.0 , 2.0, 3.0}, { 4.0, 5.0, 6.0}, {7.0, 8.0, 9.0} };
    double B[3][3] ={ {3.0 ,-2.0, 1.0}, {-2.0, 2.0, 4.0}, {3.0,-5.0, 1.0} };   

    // Print input matrices
    printf("\nInput Matrix A\n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%10.5f ", A[i][j]);
        }
        printf("\n");
    }

    printf("\nInput Matrix B\n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%10.5f ", B[i][j]);
        }
        printf("\n");
    }

    // Perform matrix addition
     mat_add(m, n, A, B, C);

    // Print output matrix
    printf("\nOutput Matrix C\n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%10.5f ", C[i][j]);
        }
        printf("\n");
    }

    return 0;
}

void mat_add(int m, int n, double A[m][n], double B[m][n], double C[m][n]) {
    for (int i = 0; i < m; i++) {
// ==============================================================================
// CODE2-1-MAT_ADD.C. A C-program for implementing Pseudocode 2.1
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
// DESCRIPTION: A module to perform C=A+B matrix addition.
//
// ON ENTRY
//   m,n :: Dimension attributes of the matrices; 
//    A  :: An input matrix (mxn);
//    B  :: An input matrix (mxn).
//
// ON EXIT
//    C :: The output matrix (mxn).
//
// REVISION DATE :: 03/18/2024
// ==============================================================================        
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}