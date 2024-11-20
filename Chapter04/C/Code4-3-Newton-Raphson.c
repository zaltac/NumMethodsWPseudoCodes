#include <math.h>
#include <stdio.h>

// ==========================================================================          
// User-defined function providing f(x), which should be cast as func(x)=0.
// ========================================================================== 
double func(double x) {
    return 4.0 + x * x * (8.0 - x * x);
}

// ==========================================================================          
// User-defined function providing f'(x)=funcp(x).
// ========================================================================== 
double funcp(double x) {
    return 4.0 * x * (4.0 - x * x);
}

// ==================================================================================
// CODE4.3-NEWTON_RAPHSON.C. A C module implementing Pseudocode 4.3.                              
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
// DESCRIPTION: A C module to compute a root of a nonlinear equation using the                   
//   Newton-Raphson method.                                                                    
//                                                                                             
// ON ENTRY                                                                                    
//   root  :: Initial guess for the root;                                                      
//   maxit :: Maximum number of iterations permitted;                                          
//   eps   :: Convergence tolerance.                                                           
//                                                                                             
// ON EXIT                                                                                     
//   iter  :: Number of iterations realized;                                                   
//   root  :: Computed approximation for the root.                                             
//                                                                                             
// USES                                                                                        
//   ABS   :: Built-in Intrinsic function returning the absolute value of a real value;        
//                                                                                             
// ALSO REQUIRED                                                                               
//   FUNC  :: User-defined external function providing the nonlinear equation, f(x).           
//   FUNCP :: User-defined external function providing the first derivative                    
//            of the nonlinear equation, f'(x).                                                
//                                                                                             
// REVISION DATE :: 11/20/2024                                                                 
// ==================================================================================
void newton_raphson(double* root, int maxit, double eps, int* iter) {
    double small, fn, fpn, aerr, rate, x0, xn, del, del0;

    small = 1.0e-15;

    printf("\n");

    del0 = 1.0;
    x0 = *root;
    *iter = 0;

    do {
        fn = func(x0);
        fpn = funcp(x0);
        del = -fn / fpn;
        aerr = fabs(del);
        rate = aerr / (del0 * del0);
        printf("%d %f %f %f %f %f\n", *iter, x0, fn, fpn, aerr, rate);
        xn = x0 + del;
        x0 = xn;
        del0 = fabs(del);
        (*iter)++;
    } while ((fabs(fn) > eps && aerr > eps) && (*iter <= maxit));

    *root = xn;

    if (*iter >= maxit) {
        printf("** Max iteration number reached=\n");
        printf("** Estimated root has NOT converged, del, f(x)= %f %f\n", fabs(del), func(xn));
    }
}

// ==============================================================================
//  Main program to test module NEWTON_RAPHSON 
// ==============================================================================
int main() {
    int maxit;
    double eps, root;
    int iter;

    maxit = 99;
    eps = 0.50e-6;
    root = 2.20;

    newton_raphson(&root, maxit, eps, &iter);

    printf("-----------------------------\n");
    printf("Root is %f converged after %d iterations\n", root, iter);
    printf("-----------------------------\n");

    return 0;
}