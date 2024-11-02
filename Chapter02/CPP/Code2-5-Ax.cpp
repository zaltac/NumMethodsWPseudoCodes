#include <iostream>
#include <vector>

std::vector<double> Ax(int n, std::vector<std::vector<double>> A, std::vector<double> x) {

// ==================================================================================
// CODE2.5-Ax.CPP. A C++ module implementing Pseudocode 2.5.                                      
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
// DESCRIPTION: A module to perform A * x = b matrix-vector multiplication.                
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Dimension attributes of input/output matrices;                                    
//    A   :: An input matrix of size nxn;                                                     
//    x   :: An input vector of length n.                                                      
//                                                                                             
// ON RETURN                                                                                     
//    b   :: The output vector of length n.                                                    
//                                                                                             
// REVISION DATE :: 03/18/2024                                                                 
// ==================================================================================
    std::vector<double> b(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * x[j];
        }
    }
    return b;
}



// ==================================================================================
// Main program to test Ax.cpp
// ==================================================================================
int main() {
    int n = 4;
    std::vector<std::vector<double>> A = {{1.0, 4.0, 3.0, -2.0},
                                         {2.0, 1.0, 2.0, 3.0},
                                         {3.0, 2.0, 1.0, 4.0},
                                         {-2.0, 3.0, 2.0, 1.0}};
    std::vector<double> x = {2.0, 5.0, 3.0, -3.0};

    std::cout << " Input Matrix A" << std::endl;
    for (const auto& row : A) {
        for (double elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n Input Vector x" << std::endl;
    for (double elem : x) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    std::vector<double> b = Ax(n, A, x);

    std::cout << "\n ------ A(n,n)X(n) product is " << std::endl;
    std::cout << "\n Output Vector b" << std::endl;

    for (double elem : b) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;


    return 0;
}