#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>


// ==============================================================================
// CODE7-2-LINEARIZE_REGRESS.CPP. A C++ program for implementing Pseudocode 7.2.
//     
// NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
// First Edition. (c) By Zekeriya ALTAÇ (2024).
// ISBN: 9781032754741 (hbk)
// ISBN: 9781032756424 (pbk)
// ISBN: 9781003474944 (ebk)
//
// DOI : 10.1201/9781003474944
// C&H/CRC PRESS, Boca Raton & London. 
//  
// This free software is complimented by the author to accompany the textbook.
// E-mail: altacz@gmail.com
//
// DESCRIPTION: A C++ module to obtain least-squares-fit to the "power model" only. Users                 
//   can likewise incorporate the other linearizable models into the module.                   
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: The number of data in the set;                                                    
//   x,y  :: Arrays of length n containing the data;                                           
//  model :: Model flag, Model = 1 corresponds to the power model, Y= a0*x^b0.                 
//                                                                                             
// ON EXIT                                                                                     
//   a0,b0:: Model parameters;                                                                 
//    E   :: Sum of the Squares of Residuals (SSR);                                            
//    S   :: Sum of the Squares of Mean Deviation (SSMD);                                      
//   r2   :: r-squared, coefficient of determination.   
//
// USES
//    exp :: Built-in Intrinsic function returning exponential of a real value, e^x.
//    pow :: Built-in intrinsic function raises a number to the power of another number.
//    log :: Built-in Intrinsic function returning the natural log of a real value.
//
// REVISION DATE :: 03/03/2025
// ==============================================================================
void LinearizeRegress(int ndata, std::vector<double>& x, std::vector<double>& y, int model, double& a0, double& b0, double& E, double& S, double& r2) {
    std::vector<double> xx(ndata), yy(ndata);
    std::vector<double> b(1, 0.0);
    std::vector<std::vector<double>> c(2, std::vector<double>(2, 0.0));
    double DD, D1, D2, yavg, yi, xk;

    for (int k = 0; k < ndata; k++) {
        if (model == 1) {
            xx[k] = std::log(x[k]);
            yy[k] = std::log(y[k]);
        } else if ((model < 1) || (model > 1)) {
            std::cerr << "Undefined model..." << std::endl;
            return;
        }
    }

    for (int i = 0; i <2; i++) {
        for (int j = 0; j <2; j++) {
            int p = i + j;
            for (int k = 0; k < ndata; k++) {
                xk = (p != 0) ? std::pow(xx[k], p) : 1.0;
                c[i][j] += xk;
            }
        }
        for (int k = 0; k < ndata; k++) {
            int p = i;
            xk = (p != 0) ? std::pow(xx[k], p) : 1.0;
            b[i] += xk * yy[k];
        }
    }
 
    DD = c[0][0] * c[1][1] - c[1][0] * c[0][1];
    D1 = b[0] * c[1][1] - b[1] * c[0][1];
    D2 = b[1] * c[0][0] - b[0] * c[1][0];

    a0 = D1 / DD;
    b0 = D2 / DD;

    yavg = 0.0;
    for (int k = 0; k < ndata; k++) {
        yavg += y[k];
    }
    yavg /= static_cast<double>(ndata);

    S = 0.0;
    E = 0.0;

    if (model <= 2) {
        a0 = std::exp(a0);
    } else {
        std::cerr << "Undefined model..." << std::endl;
        return;
    }

    for (int k = 0; k < ndata; k++) {
        S += std::pow(y[k] - yavg, 2);
        yi = a0 * std::pow(x[k], b0);
        E += std::pow(yi - y[k], 2);
    }

    r2 = 1.0 - E / S;
}


// ==============================================================================
//  The main program to test LinearizeRegress.CPP
// ==============================================================================
int main() {
    const int ndata = 6;
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<double> y = {5.3, 6.5, 6.8, 7.1, 7.8, 7.5};
    int model = 1;
    double a0, b0, E, S, r2;

    LinearizeRegress(ndata, x, y, model, a0, b0, E, S, r2);

    std::cout << "******* Best-Fit Coefficients and Parameters ********" << std::endl;
    std::cout << "  a = " << std::fixed << std::setprecision(6) << a0 << "   b = " << b0 << std::endl;
    std::cout << "  E = " << E << "   S = " << S << "   r-squared = " << r2 << std::endl;

    return 0;
}