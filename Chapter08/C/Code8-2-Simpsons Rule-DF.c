#include <stdio.h>

// ==================================================================================
// CODE8.2-Simpsons_Rule_DF.C. A C module implementing Pseudocode 8.2.                            
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
void Simpsons_Rule_DF(int n, double h, double f[], double *intg) {
    int i, m;
    double odd, even;
          
    m = n % 2;
    if (m != 0) {
        printf("Number of panels is not EVEN\n");
        return;
    }
 
    odd = 0.0;
    for (i = 1; i < n; i += 2) {
        odd += f[i];
    }

    even = 0.0;
    for (i = 2; i < n - 1; i += 2) {
        even += f[i];
    }

    *intg = f[0] + f[n] + 4.0 * odd + 2.0 * even;
    *intg = *intg * h / 3.0;
}



// ==============================================================================
//  The main program to test Simpsons_Rule_DF.C
// ==============================================================================
int main() {
    int n, i;
    double f[100], a, b, x, h, intg;

    printf("Enter Number of Panels: ");
    scanf("%d", &n);
    
    // Construct a discrete function using f(x)=x^4
    a = 0.0;
    b = 1.0;
    h = (b - a) / (double)n;
    x = a;
    
    for (i = 0; i <= n; i++) {  // Discrete dataset generation from y=FUNC(x)
        f[i] = x * x * x * x;  
        x = x + h;
    }
    
    Simpsons_Rule_DF(n, h, f, &intg);
    
    printf("%d panel Simpson's 1/3 Rule = %lf\n", n, intg);
    printf("\n");
    
    return 0;
}