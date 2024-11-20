#include <stdio.h>
#include <math.h>

double func(double x) {
// ==========================================================================          
// User-defined function providing f(x), which should be cast as func(x)=0.
// ==========================================================================  
    return x * x + 0.025 * x - 4.0;
}

// ==================================================================================
// CODE4.1-BISECTION.C. A C module implementing Pseudocode 4.1.                                   
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
// DESCRIPTION: A C module to find a root of a nonlinear equation in [a,b]                       
//   using the Bisection method.                                                               
//                                                                                             
// ON ENTRY                                                                                    
//  [a,b] :: Initial search interval (it must bracket one root);                               
//  maxit :: Maximum number of iterations permitted;                                           
//  eps   :: Convergence tolerance.                                                            
//                                                                                             
// ON EXIT                                                                                     
//  halves:: Number of halves realized;                                                        
//  root  :: Computed approximation for the root.                                              
//                                                                                             
// USES                                                                                        
//  ABS   :: Built-in Intrinsic function returning the absolute value of a real value;         
//                                                                                             
// ALSO REQUIRED                                                                               
//  FUNC  :: User-defined external function providing the nonlinear equation.                  
//                                                                                             
// REVISION DATE :: 11/20/2024                                                                 
// ==================================================================================
int Bisection(double a, double b, int maxit, double eps, double *root, int *halves) {
    double fa, fb, xm, fm, interval;
    int p = 0;
    
    fa = func(a);
    fb = func(b);
    interval = b - a;
    
    printf("p\ta\t\tb\t\tf(a)\t\tf(b)\t\txm\t\tf(xm)\t\tinterval\n");
    printf("-------------------------------------------------------------------\n");
    
    while (1) {
        p++;
        xm = 0.5 * (a + b);
        fm = func(xm);
        printf("%d\t%1.7g\t%1.7g\t%1.4e\t%1.4e\t%1.7g\t%1.4e\t%1.4e\n", p, a, b, fa, fb, xm, fm, interval);
        
        if (fa * fm > 0.0) {
            a = xm;
            fa = fm;
        } else { // case of fa*fm<0.0
            b = xm;
            fb = fm;
        }
        
        interval = 0.5 * interval;
        
        if ((fabs(fm) < eps) && (interval < eps)) {
            break;
        }
        
        if (p == maxit) {
            break;
        }
    }
    
    *root = xm;
    *halves = p;
    
    if (p == maxit) {
        printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        printf("Max iteration number reached=%d\n", maxit);
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }
    
    return 0;
}


// ==================================================================================
// Main program to test BISECTION.C
// ==================================================================================
int main() {
    int Halves, maxit = 99;
    double a = 0.0, b = 4.0, eps = 0.50e-4, root, fa, fb;
    
    fa = func(a);
    fb = func(b);
    
    if (fa * fb > 0) {
        printf("No root in interval (a,b). Change the interval.\n");
        return 1;
    }
    
    Bisection(a, b, maxit, eps, &root, &Halves);
    
    printf("\n======================================\n");
    printf("! Root is %14.8g after %d bisections !\n", root, Halves);
    printf("======================================\n");
    
    root = (a * fb - b * fa) / (fb - fa);
    printf("\n*** Root (with linear interpolation)   = %14.9g\n", root);
    printf("*** Clossness to the root, Abs[f(root)]= %10.4g\n", fabs(func(root)));
    
    return 0;
}