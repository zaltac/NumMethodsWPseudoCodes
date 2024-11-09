#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

// ==================================================================================
// CODE3.1-ENORM.CPP. A C++ module implementing ENORM module of Pseudocode 3.1.                                      
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
double ENorm(int n, vector<double> x) {
    double delta;
    delta=0.0;
    for (int i = 0; i <= n; i++) {
        delta += x[i] * x[i];
    }
    delta = sqrt(delta); 
    return delta;
}


// ==================================================================================
// CODE3.1-JACOBI.CPP. A C++ module implementing JACOBI_DRV module of Pseudocode 3.1.                                      
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
void jacobi_drv(int n, vector<vector<double>>& a, vector<double>& b, vector<double>& xo, vector<double>& x, double & delta) {
    vector<double> d(n, 0.0);  
    for (int i = 0; i < n; i++) {
        double sums = 0.0;
        for (int k = 0; k < n; k++) {
            if (i != k) {
                sums += a[i][k] * xo[k];
            }
        }
        x[i] = (b[i] - sums) / a[i][i];
        d[i] = x[i] - xo[i];
    }
    for (int k = 0; k < n; k++) {
        d[k] = x[k] - xo[k];
    }
    delta = ENorm(n,d);
}


// ==================================================================================
// CODE3.1-JACOBI.CPP. A C++ module implementing Pseudocode 3.1.                                  
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
void jacobi(int n, double eps, vector<vector<double>>& a, vector<double>& b, vector<double>& xo, int maxit, int iprint) {
    double del1, del0 = 1.0;
    int p = 0;
    vector<double> x(xo);
    while (del0 > eps || p==maxit) {
        p++;
        jacobi_drv(n, a, b, xo, x,del1);
        if (iprint == 1) {
            std::cout << "iter=" << p << "  Error= " << std::setw(10) << std::setprecision(7) << del1 << "  Ratio= " << del1 / del0 << std::endl;
        }
        xo = x;
        del0 = del1;
    }
    if (p == maxit) {
        cout << "\nJacobi method failed to converge after " << maxit << " iterations\nwithin the specified EPS tolerance." << endl;
    }
    cout << "\n------- Solution ---------------\n";
    for (int i = 0; i < n; i++) {
        cout << " x(" << i + 1 << ") = " << x[i] << endl;
    }
    cout << "\n-----------------------" << endl;
    cout << "Total no of iterations = " << p << " Maximum Error = " << del0 << endl;
}

// ==============================================================================
//  The main program to test JACOBI.CPP
// ==============================================================================
int main() {
    int n = 10;
    vector<vector<double>> a(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);
    vector<double> xo(n, 0.0);
    int iprint = 1;

    a ={ { 2.,-1., 0., 0., 0., 0., 0., 0., 0., 0.},
         {-1., 2.,-1., 0., 0., 0., 0., 0., 0., 0.},
         { 0.,-1., 2.,-1., 0., 0., 0., 0., 0., 0.},    
         { 0., 0.,-1., 2.,-1., 0., 0., 0., 0., 0.},      
         { 0., 0., 0.,-1., 2.,-1., 0., 0., 0., 0.},   
         { 0., 0., 0., 0.,-1., 2.,-1., 0., 0., 0.},
         { 0., 0., 0., 0., 0.,-1., 2.,-1., 0., 0.},    
         { 0., 0., 0., 0., 0., 0.,-1., 2.,-1., 0.},      
         { 0., 0., 0., 0., 0., 0., 0.,-1., 2.,-1.},
         { 0., 0., 0., 0., 0., 0., 0., 0.,-1., 2.} };        
    b =  { 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.2};        

    cout << "* Coefficient Matrix *" << " * RHS *" << " * Initial Guess *" << endl;

    for (int i = 0; i < n; i++) {
        std::cout << std::fixed << std::setprecision(1);
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(6) << a[i][j] << " ";
        }
        std::cout << std::setw(12) << b[i] << " " << std::setw(12) << xo[i] << std::endl;
    }    

    double eps = 1.0e-5;
    int maxit = 999;
 
    jacobi(n, eps, a, b, xo, maxit, iprint);

    return 0;
}