#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

// ==================================================================================
// CODE2.6-INV_MAT.CPP. A C++ module implementing Pseudocode 2.6.                                 
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
// DESCRIPTION: A module to find the inverse of a square matrix (with no pivoting).        
//                                                                                             
// ON ENTRY                                                                                    
//    n  :: Dimension attribute of input matrix A;                                             
//    A  :: An input matrix (nxn).                                                             
//                                                                                             
// ON RETURN                                                                                    
//    B  :: Inverse of A (nxn).                                                                
//                                                                                             
// REVISION DATE :: 03/18/2024                                                                 
// ==================================================================================
std::vector<std::vector<double>> inv_mat(std::vector<std::vector<double>> a) {
    int n = a.size();
    std::vector<std::vector<double>> b(n, std::vector<double>(n, 0.0));

    // Initialize the inverse matrix as the identity matrix
    for (int i = 0; i < n; i++) {
        b[i][i] = 1.0;
    }

    for (int j = 0; j < n; j++) {
        double p = 1.0 / a[j][j];
        for (int k = 0; k < n; k++) {
            b[j][k] *= p;
            a[j][k] *= p;
        }

        for (int i = 0; i < n; i++) {
            if (i != j) {
                double s = a[i][j];
                for (int k = 0; k < n; k++) {
                    b[i][k] -= s * b[j][k];
                    a[i][k] -= s * a[j][k];
                }
            }
        }
    }

    return b;
}


// ==================================================================================
// Main program to test Inv_Mat.cpp
// ==================================================================================
int main() {
    int n = 5;
    std::vector<std::vector<double>> a = {
        {1, 2, 1, 2, 3},
        {11., -1, 1, 4, 1},
        {4, -1, 1, 1, -1},
        {-3, 1, -8, -1, 5},
        {-1, 1, 1, 1, 1}
    };
    std::vector<std::vector<double>> b = a;

    cout << "\n ********* Input Matrix *********" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(5) << setw(10) << a[i][j] << " ";
        }
        cout << endl;
    }

    std::vector<std::vector<double>> ainv = inv_mat(a);

    cout << "\n -------- matrix A(-1) --------" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(5) << setw(10) << ainv[i][j] << " ";
        }
        cout << endl;
    }
    cout << "-------------------------------" << endl;

    std::vector<std::vector<double>> c(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                c[i][j] += b[i][k] * ainv[k][j];
            }
        }
    }

    cout << "\n -------- matrix A*A(-1) -------" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(5) << setw(10) << c[i][j] << " ";
        }
        cout << endl;
    }
    cout << "-------------------------------" << endl;

    return 0;
}
