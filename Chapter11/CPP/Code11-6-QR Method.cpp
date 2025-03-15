#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

void Basic_QR(int n, std::vector<double>& d, std::vector<double>& e, double eps, int maxit, std::vector<std::vector<double>> & V) {
// ==================================================================================
// CODE11.6-Basic_QR.cpp. A C++ module implementing Pseudocode 11.6.                              
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
//  DESCRIPTION: A C++ module implementing the QR Factorization algorithm to a symmetric       
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
//    sqrt :: Built-in CMATH library function returning the square root of a real value;           
//                                                                                             
//  REVISION DATE :: 03/15/2025                                                                
// ==================================================================================
    std::vector<double> c(n), s(n);
    for (int i = 0; i < n; i++) {
        V[i][i] = 1.0;
    }
    
    double err = 1.0;
    int p = 0;
    std::cout << " iter     Error" << std::endl;
    while (err > eps && p < maxit) {
        double t = e[0];
        for (int k = 0; k < n - 1; k++) {
            double rho = std::sqrt(d[k] * d[k] + t * t);
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
        err = std::sqrt(err);
        p++;
        std::cout << std::setw(3) << p << "    " << std::scientific << std::setprecision(4) << err << std::endl;
    }
}

int main() {
// ==============================================================================
//  The main program to test the module Basic_QR.CPP
// ==============================================================================
    int n = 4;
    std::vector<double> d = {3.0, 8.0, 6.0, 9.0};
    std::vector<double> e = {4.0, 2.0, 1.0, 0.0};
    std::vector<std::vector<double>> V(n, std::vector<double>(n, 0.0));
    
    int maxit = 500;
    double eps = 1.0e-3;
    
    Basic_QR(n, d, e, eps, maxit,V);
    
    std::cout << "\n *** Eigenvalues ***" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << " x(" << i + 1 << " )=" << std::fixed << std::setprecision(7) << d[i] << std::endl;
    }
    
    std::cout << "\n *** Eigenvectors ***" << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::fixed << std::setprecision(7) << V[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
