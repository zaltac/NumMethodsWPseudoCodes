#include <cmath>
#include <iostream>

void Mantissa_Exp(double fl, double& m, int& e) {
// ==============================================================================
// CODE1.6-MantissaExp.cpp. A C++ module implementing Pseudocode 1.6.                 
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
// DESCRIPTION: A module to determine a floating-point number's                      
//    mantissa and exponent, i.e., fl=mx10^e                                         
//                                                                                   
// ON ENTRY                                                                          
//    fl  :: A floating point number.                                                
//                                                                                   
// ON EXIT                                                                           
//    m   :: Mantissa;                                                               
//    e   :: Exponent.                                                               
//                                                                                   
// USES                                                                              
//   fabs  :: Built-in Intrinsic function returning the absolute value of a real value
//   floor :: Built-in Intrinsic function returning the greatest integer less than    
//     or equal to a real value;                                                     
//   log10 :: Built-in Intrinsic function returning the base 10 logarithm of a real value
//   pow   :: Built-in Intrinsic function returning the power of a real value.
//                                                                                   
// REVISION DATE :: 03/22/2024    
// ==============================================================================
    if (std::fabs(fl) > 0.0) {
        e = std::floor(std::log10(fl)) + 1;
        m = fl * std::pow(10.0, -e);
    } else {
        e = 0.0;
        m = 0.0;
    }
}


// ==============================================================================
//  The main program to test MANTISSA_EXP
// ==============================================================================
int main() {
    double fl, m;
    int e;

    std::cout << "Enter a real number " << std::endl;
    std::cin >> fl;

    Mantissa_Exp(fl, m, e);

    std::cout << "Mantissa=" << m << std::endl;
    std::cout << "Exponent=" << e << std::endl;

    return 0;
}