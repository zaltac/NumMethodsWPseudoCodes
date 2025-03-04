#include <iostream>
#include <vector>
#include <iomanip>


// ==================================================================================
// CODE8.2-Simpsons_Rule_DF.cpp. A C++ module implementing Pseudocode 8.2.                        
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
// DESCRIPTION: A module to estimate the integral of a discrete function f on [a,b]                      
//  using the Simpson's 1/3 rule.                                                              
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Number of panels (must be even!..);                                               
//    h   :: Interval size (uniform spacing, x_{i+1}-x_i);                                     
//    f   :: Array of length (n+1) containing the ordinates, f_0, f_1, ..., f_n.               
//                                                                                             
// ON EXIT                                                                                     
//  intg  :: Numerical estimate for the integral.                                                                        
//                                                                                             
// REVISION DATE :: 03/04/2025                                                                 
// ==================================================================================
void Simpsons_Rule_DF(int n, double h, const std::vector<double>& f, double& intg) {
    int m = n % 2;
    if (m != 0) {
        std::cout << "Number of panels is not EVEN" << std::endl;
        exit(1);
    }
    
    double odd = 0.0;
    for (int i = 1; i < n; i += 2) {
        odd += f[i];
    }
    
    double even = 0.0;
    for (int i = 2; i < n - 1; i += 2) {
        even += f[i];
    }
    
    intg = f[0] + f[n] + 4.0 * odd + 2.0 * even;
    intg = intg * h / 3.0;
}


! ==============================================================================
!  The main program to test Simpsons_Rule_DF.CPP
! ==============================================================================
int main() {
    int n;
    std::cout << "Enter Number of Panels: ";
    std::cin >> n;

    std::vector<double> f(n + 1);
    double a = 0.0;
    double b = 1.0;
    double h = (b - a) / static_cast<double>(n);
    double x = a;

    for (int i = 0; i <= n; ++i) {
        f[i] = x * x * x * x;  
        x += h;
    }

    double intg;
    Simpsons_Rule_DF(n, h, f, intg);

    std::cout << n << " panel Simpson's 1/3 Rule = " << std::setprecision(10) << intg << std::endl;
    std::cout << " " << std::endl;

    return 0;
}
