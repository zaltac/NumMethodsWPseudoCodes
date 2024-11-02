#include <stdio.h>

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
// DESCRIPTION: A module to solve a tridiagonal system of linear equations                 
//   using Thomas algorithm.                                                                   
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
void tridiagonal(int s1, int sn, double B[], double D[], double A[], double C[], double X[]) {
    for (int i = s1 + 1; i <= sn; i++) {
        double ratio = B[i] / D[i - 1];
        D[i] = D[i] - ratio * A[i - 1];
        C[i] = C[i] - ratio * C[i - 1];
    }

    X[sn - 1] = C[sn - 1] / D[sn - 1];
    for (int i = sn - 2; i >= s1; i--) {
        X[i] = (C[i] - A[i] * X[i + 1]) / D[i];
    }
}


// ==================================================================================
// Main program to test the tridiagonal.c
// ==================================================================================
int main() {
    const int n = 9;
    double b[] = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double d[] = {-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0};
    double a[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double c[] = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -14.0};
    double x[n];

    int s1 = 0;
    int sn = n;
    tridiagonal(s1, sn, b, d, a, c, x);

    printf("\n*** Solution ***\n\n");
    for (int i = 0; i < n; i++) {
        printf("x(%d)=%f\n", i + 1, x[i]);
    }
    printf("\n");

    return 0;
}