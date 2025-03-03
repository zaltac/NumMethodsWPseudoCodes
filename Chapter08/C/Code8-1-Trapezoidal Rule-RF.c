#include <stdio.h>
#include <math.h>


// ==============================================================================
// DESCRIPTION: User-defined function providing y=f(x) to be integrated. 
//
// ARGUMENTS:
//      x   :: a real input value.
// ==============================================================================
double FX(double x) {
    // User-defined function providing y=f(x)
    return pow(x, 4);
}


// ==============================================================================
// DESCRIPTION: User-defined function providing first derivative, f'(x), explicitly. 
//
// ARGUMENTS:
//      x   :: a real input value.
// ==============================================================================
double FU(double x) {
    // User-defined function providing f'(x)
    return 4.0 * pow(x, 3);
}


// ==================================================================================
// CODE8.1-Trapezoidal_Rule_RF.C. A C-module implementing Pseudocode 8.1.                         
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
// DESCRIPTION: A subroutine to estimate the integral of y=f(x) on [a,b]                       
//  using the Trapezoidal rule with/without end correction.                                    
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Number of panels (i.e., n+1 integration points);                                  
//  [a,b] :: Integration interval.                                                             
//                                                                                             
// ON EXIT                                                                                     
//  intg   :: Integral estimate using the ordinary Trapezoidal rule;                           
//  intgc  :: Integral estimate using the Trapezoidal rule with the end-point correction.      
//                                                                                             
// ALSO REQUIRED                                                                               
//    FX   :: User-defined external function providing the function, f(x);                     
//    FU   :: User-defined external function providing the first derivative, f'(x).            
//                                                                                             
// USES                                                                                        
//   double :: A built-in intrinsic function that converts an integer to a double float value. 
//                                                                                             
// REVISION DATE :: 03/03/2025                                                                 
// ==================================================================================
void Trapezoidal_Rule_RF(int n, double a, double b, double* intg, double* intgc) {
    // A subroutine to estimate the integral of y=f(x) on [a,b]
    // using the Trapezoidal rule with/without end correction.
    
    double h = (b - a) / (double)n;
    *intg = 0.50 * (FX(a) + FX(b));
    double xi = a;
    for (int i = 1; i < n; i++) {
        xi += h;
        *intg += FX(xi);
    }
    *intg *= h;
    
    double corr = -h * h * (FU(b) - FU(a)) / 12.0;
    *intgc = *intg + corr;
}


// ==============================================================================
//  The main program to test the module Trapezoidal_Rule_RF.C
// ==============================================================================
int main() {
    double a = 0.0, b = 1.0, intg, intgc;
    int n;
    
    printf("Enter Number of Panels: ");
    scanf("%d", &n);
    
    Trapezoidal_Rule_RF(n, a, b, &intg, &intgc);
    
    printf("=== Standard Trapezodial Rule ===\n");
    printf("%d, %f\n", n, intg);
    printf("\n");
    printf("=== Trapezodial Rule with End-Correction ===\n");
    printf("%d, %f\n", n, intgc);
    printf("\n");
    
    return 0;
}