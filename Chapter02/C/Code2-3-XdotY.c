#include <stdio.h>

double XdotY(int n, double* x, double* y);
// ==============================================================================
//  The main program to test FUNCTION XdotY
// ==============================================================================
int main() {
    int n = 5;
    double x[] = {1.0, 3.0, 2.0, -1.0, -1.0};
    double y[] = {-1.0, 1.0, 0.5, 1.0, 1.0};

    printf("\n Vector x ");
    for (int i = 0; i < n; i++) {
        printf("%f ", x[i]);
    }
    printf("\n Vector y ");
    for (int i = 0; i < n; i++) {
        printf("%f ", y[i]);
    }

    double dotp = XdotY(n, x, y);
    printf("\n Enorm is %f\n", dotp);

    return 0;
}

double XdotY(int n, double* x, double* y) {
// ==============================================================================
// CODE2-3-XdotY.C. A C-program for implementing Pseudocode 2.3
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
// DESCRIPTION: A function to compute the dot product of two vectors, x and y.
//
// ARGUMENTS
//    n   :: Dimension attribute of the input vectors; 
//   x, y :: The input vectors of length n.
//
// REVISION DATE :: 03/18/2024
// ==============================================================================
    double sums = 0.0;
    for (int i = 0; i < n; i++) {
        sums += x[i] * y[i];
    }
    return sums;
}