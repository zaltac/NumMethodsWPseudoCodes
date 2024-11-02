#include <iostream>
#include <cmath>
#include <vector>

std::vector<std::vector<double>> mat_add(int n, int m, const std::vector<std::vector<double>> A, const std::vector<std::vector<double>> B) {
// ==================================================================================
// CODE2.1-MAT_ADD.cpp. A C++ module implementing Pseudocode 2.1.                                 
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
//  DESCRIPTION: A module to perform C=A+B matrix addition.                                
//                                                                                             
//  ON ENTRY                                                                                   
//    m,n :: Dimension attributes of the matrices;                                             
//     A  :: An input matrix (mxn);                                                            
//     B  :: An input matrix (mxn).                                                            
//                                                                                             
//  ON RETURN                                                                                  
//     C :: The output matrix (mxn).                                                           
//                                                                                             
//  REVISION DATE :: 03/18/2024                                                                
// ==================================================================================

    std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B[i][j];
            }
        }
    return C;
}


// ==================================================================================
// Main program to test mat_add.cpp
// ==================================================================================
int main() {
    const int n = 3, m = 3;
    std::vector<std::vector<double>> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    std::vector<std::vector<double>> B = {{3.0, -2.0, 1.0}, {-2.0, 2.0, 4.0}, {3.0, -5.0, 1.0}};
    std::vector<std::vector<double>> C(m, std::vector<double>(n));

    std::cout << "\n   Input Matrix A\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::fixed << A[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\n   Input Matrix B\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::fixed << B[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    //  Perform matrix addition
    C = mat_add(m,n,A,B);


    std::cout << "\n   Output Matrix C\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::fixed << C[i][j] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}


/*

A more general alternative code

std::vector<std::vector<double>> mat_add(const std::vector<std::vector<double>> A, const std::vector<std::vector<double>> B) {

    std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));

    int mA = A.size();
    int nA = A[0].size();
    int mB = B.size();
    int nB = B[0].size();
    
    if ( mA==mB && nA==nB) {
        for (int i = 0; i < mA; i++) {
            for (int j = 0; j < nA; j++) {
                C[i][j] = A[i][j] + B[i][j];
                }
            }
        }
    else {
        std::cout << "\n   Matrices A and B are not the same size!\n";
        exit;}
    return C;
}

*/