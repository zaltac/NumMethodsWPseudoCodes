#include <cmath>
#include <iostream>

double EXPA(double x, double eps) {
// ==================================================================================
// CODE1.8-ExpA.cpp. A C++ module implementing Pseudocode 1.8.                                    
//
// NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
// First Edition. (c) By Zekeriya ALTAÇ (2024).
// ISBN: 978-1-032-75474-1 (hbk)
// ISBN: 978-1-032-75642-4 (pbk)
// ISBN: 978-1-003-47494-4 (ebk)
// 
// DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
// 
// This free software is complimented by the author to accompany the textbook.
// E-mail: altacz@gmail.com.
//                                                                                             
// DESCRIPTION: A function to compute e^x adaptively using the MacLaurin series                
//    within a user-defined tolerance.                                                         
//                                                                                             
// ARGUMENTS                                                                                   
//    x   :: A real input (exponent) value;                                                    
//   eps  :: A user-defined convergence tolerance.                                             
//                                                                                             
// USES                                                                                        
//  double:: A built-in intrinsic function that converts an integer argument to a double precision value.  
//  fabs  :: Built-in Intrinsic function returning the absolute value of a real value.            
//                                                                                             
// REVISION DATE :: 04/11/2024                                                                                                                                                       
// ==================================================================================
    double sum = 1.0, term = 1.0;
    int k = 0;
    while (fabs(term) > eps * fabs(sum)) {
        k++;
        term *= x / static_cast<double>(k);
        sum += term;
    }
    return sum;
}


// ==============================================================================
//  The main program to test function EXPA
// ==============================================================================
int main() {
    double x, eps = 1.0e-4;
    std::cout << "Enter x" << std::endl;
    std::cin >> x;
    std::cout << "e^" << x << " = " << EXPA(x, eps) << std::endl;
    return 0;
}