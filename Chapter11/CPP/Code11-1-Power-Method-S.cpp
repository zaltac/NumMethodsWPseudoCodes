#include <iostream>
#include <cmath>
#include <vector>

double MAX_SIZE(int n, const std::vector<double>& X) {
    // ==================================================================================
    // CODE11.1-MAX_SIZE.CPP of a module in Pseudocode 11.1.                            
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
    //  DESCRIPTION: A function module to find the largest element (in absolute value) of an array.
    //                                                                                             
    //  ARGUMENTS                                                                                  
    //     n   :: Length of the array;                                                          
    //     x   :: Array of length n.                                                            
    //                                                                                             
    //  USES                                                                                       
    //    abs  :: Built-in CMATH library function returning the absolute value of a real value;       
    //                                                                                             
    //  REVISION DATE :: 03/10/2025                                                                
    // ==================================================================================
    double xmax = X[0];
    for (int i = 1; i < n; i++) {
        if (abs(X[i]) > abs(xmax)) {
            xmax = X[i];
        }
    }
    return xmax;
}

double ENORM(int n, const std::vector<double>& c) {
    // ==================================================================================
    // CODE3.1-ENORM.CPP of a module in Pseudocode 3.1.   
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
    // DESCRIPTION: A C++ function module to compute Euclidean (L2) norm of a vector.                  
    //                                                                                             
    // ARGUMENTS                                                                                   
    //     n  :: The length of an input vector;                                                    
    //     x  :: A vector (array) of length n.                                                     
    //                                                                                             
    // USES                                                                                        
    //   sqrt :: Built-in CMATH library function returning the square root of a real value.            
    //                                                                                             
    // REVISION DATE :: 11/09/2024                                                                  
    // ==================================================================================
    double ENORM = 0.0;
    for (int k = 0; k < n; k++) {
        ENORM += c[k] * c[k];
    }
    return sqrt(ENORM);
}

void AX(int n, const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& b) {
    // ==================================================================================
    // CODE2.5-AX.CPP of a module in Pseudocode 2.5.                                        
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
    // DESCRIPTION: A C++ module to perform A * x = b matrix-vector multiplication.                
    //                                                                                             
    // ON ENTRY                                                                                    
    //    n   :: Dimension attributes of input/output matrices;                                    
    //    A   :: An input matrix of size nxn;                                                     
    //    x   :: An input vector of length n.                                                      
    //                                                                                             
    // ON EXIT                                                                                     
    //    b   :: The output vector of length n.                                                    
    //                                                                                             
    // REVISION DATE :: 03/18/2024                                                                 
    // ==================================================================================
    for (int i = 0; i < n; i++) {
        b[i] = 0.0;
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * x[j];
        }
    }
}

void POWER_METHOD_S(int n, const std::vector<std::vector<double>>& A, std::vector<double>& X, double* eigen, double* Error, double eps, int maxit) {
    // ==================================================================================
    // CODE11.1-POWER_METHOD_S.cpp. A C++ module implementing Pseudocode 11.1.                        
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
    //  DESCRIPTION: A C++ module to find the dominant eigenvalue (largest in absolute value)          
    //     using the Power Method with scaling technique.                                          
    //                                                                                             
    //  ON ENTRY                                                                                   
    //     n     :: Size of the matrix;                                                            
    //     A     :: A real square matrix (nxn);                                                    
    //     x     :: Array of length n containing the initial guess for the eigenvector;            
    //    lambda :: An initial guess for the dominant eigenvalue;                                  
    //     eps   :: Convergence tolerance;                                                         
    //    maxit  :: Maximum number of iterations permitted.                                        
    //                                                                                             
    //  ON EXIT                                                                                    
    //    lambda :: Estimated dominant (largest in absolute value) eigenvalue;                     
    //     x     :: Array of length n containing the estimated eigenvector;                        
    //    error  :: Error, max of both L2-norm of displacement vector and relative error           
    //               for eigenvalue.                                                               
    //                                                                                             
    //  USES                                                                                       
    //    abs     :: Built-in CMATH library function returning the absolute value of a real value;
    //    AX      :: A module for computing Ax matrix-vector product (see Pseudocode 2.5);          
    //    MAX_SIZE:: A function module providing largest (in absolute) value of a vector.          
    //                                                                                             
    //  REVISION DATE :: 03/10/2025                                                                
    // ==================================================================================
    std::vector<double> d(n), Xn(n);
    double Err1, Err2, eigeno;
    int p = 0;
    
    eigeno = *eigen;

    *Error = 0.0;
    for (int i = 0; i < n; i++) {
        if (fabs(X[i]) > fabs(*Error)) {
            *Error = X[i];
        }
    }

    while (*Error > eps && p < maxit) {
        p++;
        AX(n, A, X, Xn);

            // Find max of |Xn|
        *eigen = MAX_SIZE(n, Xn);
            // Normalize Xn
        for (int i = 0; i < n; i++) {
            Xn[i] /= *eigen;
        }

        for (int i = 0; i < n; i++) {
            d[i] = fabs(Xn[i] - X[i]);
        }
        Err1 = ENORM(n, d);
        Err2 = fabs(1.0 - eigeno / *eigen);
        *Error = Err1 > Err2 ? Err1 : Err2;

        std::cout << "iter=" << p << " eigenval=" << *eigen << " Abs error=" << *Error << std::endl;
        std::cout << "Eigenvector: ";
        for (int i = 0; i < n; i++) {
            std::cout << Xn[i] << " ";
        }
        std::cout << std::endl << std::endl; 

        eigeno = *eigen;
        for (int i = 0; i < n; i++) {
            X[i] = Xn[i];
        }
    }
}



int main() {
    // ==============================================================================
    //  The main program to test the module POWER_METHOD_S.CPP
    // ==============================================================================
    int n = 4;
    std::vector<std::vector<double>> A = {{  8.0,   9.0,  10.0,   9.0},
                                          {-13.0, -12.0, -12.0, -11.0},
                                          {-18.0,  -9.0, -20.0,  -9.0},
                                          { 11.0,   1.0,  10.0,   0.0}};
    std::vector<double> X = {1.0, 0.0, 0.0, 0.0};
    double eigen, Aerr, eps = 1.0e-3;
    int maxit = 300;

    // Apply Power Method
    POWER_METHOD_S(n, A, X, &eigen, &Aerr, eps, maxit);

    if (Aerr < eps) {
        std::cout << "===============================\n";
        std::cout << "Largest eigenvalue=" << eigen << std::endl;
        std::cout << "Eigenvector       =";
        for (int i = 0; i < n; i++) {
            std::cout << X[i] << " ";
        }
        std::cout << "\n===============================\n";
    } else {
        std::cout << "Max number of iterations is reached.\n";
    }

    return 0;
}