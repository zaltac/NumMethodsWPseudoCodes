#include <stdio.h>

double macheps(double x); // Declare macheps

// ==============================================================================
//  The main program to MACHEPS
// ==============================================================================
int main() {
    double x, y;
    printf("Enter a real value in 0<x<=1\n");
    scanf("%lf", &x);
    
    // Print result
    y = macheps(x);
    
    printf("Machine Epsilon is = %e \n",y);
    return 0;
}

double macheps(double x) {
// ==================================================================================
// CODE1.5-MachEps.c. A C module implementing Pseudocode 1.5 in                     
//
// NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
// First Edition. (c) By Zekeriya ALTA� (2024).
// ISBN: 978-1-032-75474-1 (hbk)
// ISBN: 978-1-032-75642-4 (pbk)
// ISBN: 978-1-003-47494-4 (ebk)
// 
// DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
// 
// This free software is complimented by the author to accompany the textbook.
// E-mail: altacz@gmail.com.
//                                                                                   
// DESCRIPTION: A function module to calculate the machine epsilon, which            
//    provides a positive machine value that is almost negligible compared to 1.     
//                                                                                   
// INPUT ARGUMENT                                                                    
//    x   :: A real number, 0<x<=1.                                                  
//                                                                                   
// REVISION DATE :: 03/18/2024                                                       
//                                                                                   
// ==================================================================================  
    while (1.0 + x/2.0 > 1.0) {
        x = x/ 2.0;
    };
    return x;
}

