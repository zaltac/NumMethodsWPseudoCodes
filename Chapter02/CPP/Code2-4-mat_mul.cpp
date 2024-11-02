#include <iostream>
#include <vector>

std::vector<std::vector<double>> mat_mul(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
// ==================================================================================
// CODE2.4-MAT_MUL.CPP. A C++ module implementing Pseudocode 2.4.                                     
//
// NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
// First Edition. (c) By Zekeriya ALTA� (2024).
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
//  DESCRIPTION: A module to find A*B=C matrix multiplication.                             
//                                                                                             
//  ON ENTRY                                                                                   
//   m,p,n :: Dimension attributes of input/output matrices;                                   
//      A  :: An input matrix of size m×p;                                                    
//      B  :: An input matrix of size p×n.                                                    
//                                                                                             
//  ON RETURN                                                                                    
//      C  :: The output matrix of size m×n.                                                  
//                                                                                             
//  REVISION DATE :: 03/18/2024                                                                
// ==================================================================================
    int m = A.size();
    int pA = A[0].size();
    int pB = B.size();
    int n = B[0].size();

    std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));

    if (pA == pB) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < pA; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

    }
    else {
        std::cout << "Sorry, cannot multiply A and B." << std::endl;
        std::cout << "Execution halted!!!" << std::endl;
        exit(0);
    }

    return C;
}

// ==============================================================================
//  The main program to test mat_mul.cpp
// ==============================================================================
int main() {
    std::vector<std::vector<double>> A = {{1.0, -3.0, 5.0, 1.0}, 
                                         {-2.0, 4.0, 1.0, 2.0}};

    std::vector<std::vector<double>> B = {{3.0, 1.0, -2.0}, 
                                         {2.0, 0.0, 4.0}, 
                                         {-1.0, 1.0, -3.0}, 
                                         {2.0, 5.0, 3.0}};

    std::cout << "Input Matrix A \n" << std::endl;
    for (const auto& row : A) {
        for (double elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\n Input Matrix B \n" << std::endl;
    for (const auto& row : B) {
        for (double elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    std::vector<std::vector<double>> C = mat_mul(A, B);

    std::cout << "\n ------ A(m,p)xB(p,n) matrix product is ------ \n" << std::endl;
    for (const auto& row : C) {
        for (double elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\n ------ Matrix product with std::matmul ------ \n" << std::endl;
    C = mat_mul(A, B);
    for (const auto& row : C) {
        for (double elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}