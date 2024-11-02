#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

// Function to solve an upper triangular system using forward substitution
vector<double> back_substitute(int n, vector<vector<double>>& a, vector<double>& b) {
    vector<double> x(n);
    x[n - 1] = b[n - 1] / a[n - 1][n - 1];
    for (int k = n - 2; k >= 0; k--) {
        double sums = 0.0;
        for (int j = k + 1; j < n; j++) {
            sums += a[k][j] * x[j];
        }
        x[k] = (b[k] - sums) / a[k][k];
    }
    return x;
}

// ==================================================================================
// CODE2.9-LINEAR_SOLVE.cpp. A C++ module implementing Pseudocode 2.9.                            
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
// DESCRIPTION: A module to solve a system of linear equations using naive                 
//   Gauss Elimination (opt=0) or Gauss-Jordan Elimination (opt/=0) algorithm.                 
//                                                                                             
// ON ENTRY                                                                                    
//     n  :: Number of unknowns;                                                               
//     A  :: Input coefficient matrix of size nÃ—n;                                            
//     b  :: An input array of length n containing the rhs;                                    
//    opt :: Option key (=0, Naive Gauss-Elim.; /=0, Gauss-Jordan Elimn).                      
//                                                                                             
// ON RETURN                                                                                     
//     x  :: The output array of length n containing the solution.                             
//                                                                                             
// USES                                                                                        
//   fabs  :: Built-in Intrinsic function returning the absolute value of a real value;         
//   back_substitute :: a module to solve an upper-triangular system.                    
//                                                                                             
// REVISION DATE :: 03/18/2024                                                                 
// ==================================================================================
vector<double> linear_solve(int n, vector<vector<double>>& a, vector<double>& b, int opt) {
    double eps = 1e-12;

    // Forward elimination steps
    for (int j = 0; j < n; j++) {
        double ajj = a[j][j];
        if (fabs(ajj) < eps) {
            cout << "Pivot is zero at j=" << j << endl;
            cout << "Execution is halted!" << endl;
            return {}; // return empty vector
        } else {
            double val = 1.0 / ajj;
            b[j] *= val;
            for (int k = 0; k < n; k++) {
                a[j][k] *= val;
            }
            for (int i = j + 1; i < n; i++) {
                double s = a[i][j];
                a[i][j] = 0.0;
                for (int k = j + 1; k < n; k++) {
                    a[i][k] -= s * a[j][k];
                }
                b[i] -= s * b[j];
            }
        }
    }

    // Back substitution steps
    vector<double> x(n);
    if (opt == 0) {
        x = back_substitute(n, a, b);
    } else {
        for (int j = n - 1; j >= 0; j--) {
            for (int i = j - 1; i >= 0; i--) {
                b[i] -= a[i][j] * b[j];
                a[i][j] = 0.0;
            }
            x[j] = b[j];
        }
        x[0] = b[0];
    }

    return x;
}


// ==================================================================================
// Main program to test the linear_solve.cpp
// ==================================================================================
int main() {
    int n = 3;
    vector<vector<double>> a = {{1., 4., -6.}, {-1., 6., -4.}, {4., -1., -1}};
    vector<double> b = {0., 60., 0.};

    int opt = 0;
    vector<double> x = linear_solve(n, a, b, opt);

    cout << "\n Method applied here is" << endl;
    if (opt == 0) {
        cout << " Naive Gauss Elimination" << endl;
    } else {
        cout << " Gauss Jordan Elimination" << endl;
    }

    cout << "\n------ Matrix A & b ------------\n" << endl;

    for (int i = 0; i < n; i++) {
        std::cout << std::fixed << std::setprecision(3);
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(8) << a[i][j] << " ";
        }
        std::cout << std::setw(12) << b[i] << " " << std::endl;
    } 

    cout << "\n------- Solution ---------------\n" << endl;
    for (int i = 0; i < n; i++) {
        cout << " x(" << i + 1 << ") = " << x[i] << endl;
    }
 
    return 0;
}
