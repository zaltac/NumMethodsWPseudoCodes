#include <stdio.h>
#include <math.h>

#define N 4

void Basic_QR(int n, double d[], double e[], double eps, int maxit, double V[N][N]) {
// ==================================================================================
// CODE11.6-Basic_QR.c. A C module implementing Pseudocode 11.6.                              
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
//  DESCRIPTION: A C module implementing the QR Factorization algorithm to a symmetric       
//     tridiagonal matrix to find its eigenvalues and eigenvectors.                            
//                                                                                             
//  ON ENTRY                                                                                   
//     n   :: Dimension attribute of the tridiagonal matrix (nxn);                             
//     d   :: An array of length n containing the main diagonal, d(1) ... d(n);                
//     e   :: An array of length n containing the subdiagonal, e(1) ... e(n-1).                
//                                                                                             
//  ON EXIT                                                                                    
//     d   :: An array of length n containing the eigenvalues;                                 
//     V   :: A square matrix containing the eigenvector(nxn).                                 
//                                                                                             
//  USES                                                                                       
//    sqrt :: Built-in MATH library function returning the square root of a real value;           
//                                                                                             
//  REVISION DATE :: 03/15/2025                                                                
// ==================================================================================
    double c[n], s[n];
    for (int i = 0; i < n; i++) {
        V[i][i] = 1.0;
    }
    
    double err = 1.0;
    int p = 0;
    printf(" iter     Error\n");
    while (err > eps && p < maxit) {
        double t = e[0];
        for (int k = 0; k < n - 1; k++) {
            double rho = sqrt(d[k] * d[k] + t * t);
            c[k] = d[k] / rho;
            s[k] = t / rho;
            d[k] = rho;
            t = e[k];
            e[k] = t * c[k] + d[k + 1] * s[k];
            d[k + 1] = -t * s[k] + d[k + 1] * c[k];
            if (k < n - 2) {
                t = e[k + 1];
                e[k + 1] = t * c[k];
            }
            for (int i = 0; i < n; i++) {
                double q = V[i][k];
                double r = V[i][k + 1];
                V[i][k] = c[k] * q + s[k] * r;
                V[i][k + 1] = -s[k] * q + c[k] * r;
            }
        }
        
        for (int k = 0; k < n - 1; k++) {
            d[k] = d[k] * c[k] + e[k] * s[k];
            double t = d[k + 1];
            e[k] = t * s[k];
            d[k + 1] = t * c[k];
        }
        
        err = 0.0;
        for (int i = 0; i < n; i++) {
            err += e[i] * e[i];
        }
        err = sqrt(err);
        p++;
        printf("%3d    %e\n", p, err);
    }
}

int main() {
// ==============================================================================
//  The main program to test the module Basic_QR.C
// ==============================================================================
    double d[N] = { 3.0, 8.0, 6.0, 9.0};
    double e[N] = { 4.0, 2.0, 1.0, 0.0};
    double V[N][N] = {0};
    
    int maxit = 500;
    double eps = 1.0e-3;
    
    Basic_QR(N, d, e, eps, maxit, V);
    
    printf("\n *** Eigenvalues ***\n");
    for (int i = 0; i < N; i++) {
        printf(" x(%d)=%.7f\n", i + 1, d[i]);
    }
    
    printf("\n *** Eigenvectors ***\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.7f ", V[i][j]);
        }
        printf("\n");
    }
    return 0;
}