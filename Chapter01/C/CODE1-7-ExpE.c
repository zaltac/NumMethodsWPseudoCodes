#include <stdio.h>

// Declare EXPE
double EXPE(double x, int n);


// ==============================================================================
//  The main program to test function EXPE
// ==============================================================================
int main() {
    double x;
    int n;
    
    printf("Enter x, n\n");
    scanf("%lf %d", &x, &n);
    
    double result = EXPE(x, n);
    printf("e**%lf = %.6f\n", x, result);
    
    return 0;
}

double EXPE(double x, int n) {  
// ==================================================================================
// CODE1.7-ExpE.c. A C function module implementing Pseudocode 1.7.                                      
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
// DESCRIPTION: A function to compute e^x using the MacLaurin series with specified            
//    number of terms.                                                                         
//                                                                                             
// ARGUMENTS                                                                                   
//    x   :: A real input (exponent) value;                                                    
//    n   :: The number of terms of the MacLauring series to be included.                      
//                                                                                             
// USES                                                                                        
//   double:: A built-in intrinsic function that converts an integer argument to a double value.  
//                                                                                             
// REVISION DATE :: 03/18/2024                                                                                                                                                              
// ==================================================================================   
    double sum = 1.0;
    double term = 1.0;
    
    for (int k = 1; k < n; k++) {
        term *= x / (double)k;   // (double)k converts k to double precision float
        sum += term;
    }
    
    return sum;
}