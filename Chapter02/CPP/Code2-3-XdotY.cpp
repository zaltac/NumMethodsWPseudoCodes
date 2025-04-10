#include <iostream>
#include <vector>

double XdotY(int n, std::vector<double>& x, std::vector<double>& y) {
// ==================================================================================
// CODE2.3-XdotY.cpp. A C++ module implementing Pseudocode 2.3.                                   
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
//  DESCRIPTION: A function to compute the dot product of two vectors, x and y.                
//                                                                                             
//  ARGUMENTS                                                                                  
//     n   :: Dimension attribute of the input vectors;                                        
//    x, y :: The input vectors of length n.                                                   
//                                                                                             
//  REVISION DATE :: 03/18/2024                                                                
// ==================================================================================
    double sums = 0.0;
    for (int i = 0; i < n; i++) {
        sums += x[i] * y[i];
    }
    return sums;
}


// ==================================================================================
// Main program to test XdotY.cpp
// ==================================================================================
int main() {
    int n = 5;
    std::vector<double> x = {1.0, 3.0, 2.0, -1.0, -1.0};
    std::vector<double> y = {-1.0, 1.0, 0.5, 1.0, 1.0};

    std::cout << "\n Vector x ";
    for (double value : x) {
        std::cout << value << " ";
    }
    std::cout << "\n Vector y ";
    for (double value : y) {
        std::cout << value << " ";
    }

    double dotp = XdotY(n, x, y);
    std::cout << "\n Enorm is " << dotp << std::endl;

    return 0;
}