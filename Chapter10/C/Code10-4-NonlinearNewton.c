#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double ENORM(int n, double c[]) {
// ==================================================================================
// CODE3.1-ENORM.C. A C module implementing ENORM of Pseudocode 3.1.                                      
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
// DESCRIPTION: A function module to compute Euclidean (L2) norm of a vector.                  
//                                                                                             
// ARGUMENTS                                                                                   
//     n  :: The length of an input vector;                                                    
//     x  :: A vector (array) of length n.                                                     
//                                                                                             
// USES                                                                                        
//   sqrt :: Built-in intrinsic function returning the square root of a real value.         
//                                                                                             
// REVISION DATE :: 11/09/2024                                                                  
// ==================================================================================
    double norm = 0.0;
    for (int k = 0; k < n; ++k) {
        norm += c[k] * c[k];
    }
    return sqrt(norm);
}

void FUNCS(double x, double y, double yp, double* fun, double* dfdy, double* dfdp) {
   // ==============================================================================
   // DESCRIPTION: A user-defined function supplying the coefficients of the nonlinear 
   //    two-point BVP given in the form: 
   //              y'' = f(x,y,y')  on [a,b]
   //
   // ON ENTRY
   //   x    :: Independent variable (a <=x<= b);
   //   y    :: Dependent variable y=y(x).
   //
   // ON EXIT
   //   f    :: The nonlinear two-point BVP f(x,y,y') evaluated at (x,y);
   //   yp   :: First derivative of the dependent variable y';
   //   dfdy :: Partial derivative of f wrt y evaluated at (x,y), df/dy;
   //   dfdp :: Partial derivative of f wrt y' evaluated at (x,y), df/dy'.
   //
   // REVISION DATE :: 03/10/2025
   // ==============================================================================
    *fun = y - yp*yp/y;
    *dfdy = 1.0 + (yp/y)*(yp/y);
    *dfdp = -2.0*yp/y;
}

double Exact(double x) {
   // ==============================================================================
   // DESCRIPTION: A function program providing the true solution y=f(x) for 
   //    testing the example problem. 
   //
   // ARGUMENTS:
   //      x   :: A double input, independent variable.
   //
   // USES                                                                                        
   //   sqrt :: Built-in intrinsic function returning the square root of a real value;
   //   cosh :: Built-in Intrinsic function returning the hyperbolic cosine of a real value.  
   // ==============================================================================
    double term = cosh(sqrt(2.0) * (1.0 - x));
    term = sqrt(term / cosh(sqrt(2.0)));
    return term;
}

void TRIDIAGONAL(int s1, int sn, double B[], double D[], double A[], double C[], double X[]) {
   // ==================================================================================
   // CODE2.13-TRIDIAGONAL.c. A C module implementing Pseudocode 2.13.                               
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
   // DESCRIPTION: A module to solve a tridiagonal system of linear equations                 
   //   using Thomas algorithm.                                                                   
   //                                                                                             
   // ON ENTRY                                                                                    
   //    s1 :: Subscript of the first unknown (usually 1);                                        
   //    sn :: Subscript of the last unknown (usually no. of eqs, n);                             
   //     b :: Array of length n containing coefficients of below diagonal elements;              
   //     d :: Array of length n containing coefficients of diagonal elements;                    
   //     a :: Array of length n containing coefficients of above diagonal elements;              
   //     c :: Array of length n containing coefficients of rhs.                                  
   //                                                                                             
   // ON EXIT                                                                                     
   //     x :: An array of length n containing the solution.                                      
   //                                                                                             
   // REVISION DATE :: 03/18/2024                                                                 
   // ==================================================================================
    for (int i = s1; i < sn; ++i) {
        double ratio = B[i] / D[i - 1];
        D[i] -= ratio * A[i - 1];
        C[i] -= ratio * C[i - 1];
    }

    X[sn - 1] = C[sn - 1] / D[sn - 1];
    for (int i = sn - 2; i >= s1; --i) {
        X[i] = (C[i] - A[i] * X[i + 1]) / D[i];
    }
}


void NONLINEAR_NEWTON(int M, double x[], double yo[], double y[], double tol, double alpha[], double beta[], double gamma[], int maxit) {
// ==================================================================================
// CODE10.4-NONLINEAR_NEWTON.C. A C module implementing Pseudocode 10.4.                          
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
//  DESCRIPTION: A module to find approximate solution of a two-point nonlinear differential   
//    equation using the Newton's Method. The nonlinear equations is cast in the following form
//          y'' = f(x,y,y')  on [a,b]                                                          
//    subject to                                                                               
//          alpha1 * y'(a)+ beta1 * y(a) = gamma1                                              
//          alpha2 * y'(b)+ beta2 * y(b) = gamma2                                              
//                                                                                             
//  CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1. 
//                                                                                             
//  ON ENTRY                                                                                   
//     M   :: Number of (equations) grid poinds;                                               
//     x   :: Array of length M containing the abscissas of the grid points;                   
//     yo  :: Array of length M containing the initial guess for the solution;                 
//    alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                  
//           the boundary conditions as stated above;                                          
//    eps  :: Convergence tolerance;                                                           
//    maxit:: Maximum number of iterations permitted.                                          
//                                                                                             
//  ON EXIT                                                                                    
//     y   :: Array of length M containing the approximate solution.                           
//                                                                                             
//  USES                                                                                       
//    FUNCS :: A user-defined external function module providing the coefficients of the 
//            nonlinear two-point BVP;           
//    ENORM:: A function module to calculate the Euclidean vector (L2 norm) of a vector;       
//    TRIDIAGONAL :: A module to solve a tridiagonal system of equations with Thomas algorithm.
//                                                                                             
//  REVISION DATE :: 03/10/2025                                                                
// ==================================================================================
    double A[M], B[M], D[M], C[M];
    double h = x[1] - x[0];
    double hsqr = h * h;
    double bb2h = 0.5 / h;
    double hb2 = 0.5 * h;
    double c1d = 2.0 * h * beta[0] / alpha[0];
    double c1r = 2.0 * h * gamma[0] / alpha[0];
    double c2d = 2.0 * h * beta[1] / alpha[1];
    double c2r = 2.0 * h * gamma[1] / alpha[1];

    int bc_left = abs(alpha[0]);
    int bc_rigt = abs(alpha[1]);

    double aerr = 1.0;
    int p = 0;
    for (int i = 0; i < M; ++i) {
        y[i] = yo[i];
    }

    while (aerr > tol && p < maxit) {
        for (int k = 0; k < M; ++k) {
            if (k == 0) {
                if (bc_left == 0) {
                    D[k] = 1.0;
                    yo[k] = gamma[0] / beta[0];
                    A[k] = 0.0;
                    B[k] = 0.0;
                    C[k] = 0.0;
                } else {
                    double yx = (gamma[0] - beta[0] * yo[k]) / alpha[0];
                    double fun, dfdy, dfdp;
                    FUNCS(x[k], yo[k], yx, &fun, &dfdy, &dfdp);
                    A[k] = 2.0;
                    B[k] = 0.0;
                    D[k] = -2.0 + c1d - hsqr * (dfdy - dfdp * beta[0] * yo[k] / alpha[0]);
                    C[k] = -c1r + (-2.0 + c1d) * yo[k] + 2.0 * yo[k + 1] - hsqr * fun;
                }
            } else if (k == M - 1) {
                if (bc_rigt == 0) {
                    D[k] = 1.0;
                    yo[k] = gamma[1] / beta[1];
                    A[k] = 0.0;
                    B[k] = 0.0;
                    C[k] = 0.0;
                } else {
                    double yx = (gamma[1] - beta[1] * yo[k]) / alpha[1];
                    double fun, dfdy, dfdp;
                    FUNCS(x[k], yo[k], yx, &fun, &dfdy, &dfdp);
                    B[k] = 2.0;
                    D[k] = -2.0 - c2d - hsqr * (dfdy - dfdp * beta[1] * yo[k] / alpha[1]);
                    A[k] = 0.0;
                    C[k] = c2r - (2.0 + c2d) * yo[k] + 2.0 * yo[k - 1] - hsqr * fun;
                }
            } else {
                double yx = bb2h * (yo[k + 1] - yo[k - 1]);
                double yxx = yo[k + 1] - 2.0 * yo[k] + yo[k - 1];
                double fun, dfdy, dfdp;
                FUNCS(x[k], yo[k], yx, &fun, &dfdy, &dfdp);
                B[k] = 1.0 + hb2 * dfdp;
                D[k] = -2.0 - hsqr * dfdy;
                A[k] = 1.0 - hb2 * dfdp;
                C[k] = yxx - hsqr * fun;
            }
        }

        TRIDIAGONAL(1, M, B, D, A, C, C);
        aerr = ENORM(M, C);
        for (int k = 0; k < M; ++k) {
            y[k] = yo[k] - C[k];
        }
        p++;
        for (int k = 0; k < M; ++k) {
            yo[k] = y[k];
        }
    }
}



int main() {
// ==============================================================================
//  The main program to test the module NONLINEAR_NEWTON.C
// ==============================================================================
    #define nmax 201

    int n;
    printf("Enter No of intervals: ");
    scanf("%d", &n);
    if (n > nmax) {
        printf("Max grid size is %d\n", nmax);
        printf("Increase NMAX, and then try it again.\n");
        return 1;
    }

    double tol = 1e-6;
    double guess = 0.4;
    int maxit = 99;
    int neq = n + 1;
    double xa = 0.0, xb = 1.0;
    double h = (xb - xa) / (double)n;

    double x[neq], y[neq], yo[neq], alpha[2], beta[2], gamma[2];

    for (int i = 0; i < neq; ++i) {
        x[i] = xa + (double)i * h;
    }

    alpha[0] = 0.0; beta[0] = 1.0; gamma[0] = 1.0;
    alpha[1] = 1.0; beta[1] = 0.0; gamma[1] = 0.0;

    int bc_left = abs(alpha[0]);
    int bc_rigt = abs(alpha[1]);

    for (int i = 0; i < neq; ++i) {
        x[i] = xa + (double)i * h;
        yo[i] = guess;
        if (bc_left == 0 && bc_rigt == 0) {
            yo[i] = gamma[0] + (gamma[1] - gamma[0]) * (x[i] - xa) / (xb - xa);
        }
    }

    NONLINEAR_NEWTON(neq, x, yo, y, tol, alpha, beta, gamma, maxit);

    printf("\n        x        Exact     N.Approx   Abs Error\n");
    for (int i = 0; i < neq; ++i) {
        double aerr = fabs(Exact(x[i]) - y[i]);
        printf("%10.5f%12.7f%12.7f%15.5e\n", x[i], Exact(x[i]), y[i], aerr);
    }

    printf("\n*** D O N E ***\n");

    return 0;
}