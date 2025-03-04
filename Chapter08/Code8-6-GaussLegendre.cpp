#include <cmath>
#include <iostream>
#include <iomanip>


// ==================================================================================
// CODE8.6-Gauss_Legendre_Quad.cpp. A C++ module implementing Pseudocode 8.6.                     
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
// DESCRIPTION: A subroutine to generate N-point Gauss-Legendre quadrature                     
//   abscissas and weights on [-1,1].                                                          
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Number of quadrature points;                                                      
//    eps :: Tolerance, i.e., desired level of numerical accuracy.                             
//                                                                                             
// ON EXIT                                                                                     
//    x   :: Array of length N containing the abscissas;                                       
//    w   :: Array of length N containing the weights.                                         
//                                                                                             
// USES                                                                                        
//   fabs :: Built-in Intrinsic function returning the absolute value of a real value;         
//   cos  :: Built-in Intrinsic function returning trig cosine value.                          
//                                                                                             
// REVISION DATE :: 03/04/2025                                                                 
// ==================================================================================
void Gauss_Legendre_Quad(int n, double eps, double x[], double w[]) {
    int m, i, k;
    double u, del, P0, P1, P2, PP, pi = 3.1415926535897932385;
    
    m = (n + 1) / 2;
    for (i = 0; i < m; i++) {
        u = cos(pi * (4 * i + 3) / (4 * n + 2));
        del = 1.0;
        while (fabs(del) > eps) {
            P0 = 1.0;
            P1 = u;
            for (k = 2; k <= n; k++) {
                P2 = (2 * k - 1) * u * P1 - (k - 1) * P0;
                P2 /= k;
                P0 = P1;
                P1 = P2;
            }
            PP = n * (u * P2 - P0) / (u * u - 1.0);
            del = P1 / PP;
            u -= del;
        }
        x[i] = -u;
        w[i] = 2.0 / (1.0 - u * u) / (PP * PP);
        x[n - i - 1] = u;
        w[n - i - 1] = w[i];
    }
}


// ==============================================================================
//  The main program to test SUBROUTINE Gauss_Legendre_Quad.CPP
// ==============================================================================
int main() {
    int n, i;
    double x[90], w[90];
    
    std::cout << "Enter n" << std::endl;
    std::cin >> n;
    
    double eps = 1.0e-6;
    Gauss_Legendre_Quad(n, eps, x, w);
    
    std::cout << std::setw(5) << "i" << std::setw(18) << "x_i" << std::setw(18) << "w_i" << std::endl;
    std::cout << std::setw(5) << "---" << std::setw(18) << "---------------" << std::setw(18) << "---------------" << std::endl;
    for (i = 0; i < n; i++) {
        std::cout << std::setw(5) << i + 1 << std::setw(18) << std::fixed << std::setprecision(12) << x[i] << std::setw(18) << std::fixed << std::setprecision(12) << w[i] << std::endl;
    }
    
    return 0;
}