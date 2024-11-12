#include <stdio.h>
#include <math.h>


// ==================================================================================
// CODE2.5-Ax.C. A C module implementing Pseudocode 2.5.                                          
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
// DESCRIPTION: A module to perform A * x = b matrix-vector multiplication.                
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Dimension attributes of input/output matrices;                                    
//    A   :: An input matrix of size nxn;                                                     
//    x   :: An input vector of length n.                                                      
//                                                                                             
// ON RETURN                                                                                     
//    b   :: The output vector of length n.                                                    
//                                                                                             
// REVISION DATE :: 03/18/2024                                                                 
// ================================================================================== 
void AX(int n, double A[][4], double x[], double c[]) {
    for (int i = 0; i < n; i++) {
        c[i] = 0.0;
        for (int j = 0; j < n; j++) {
            c[i] += A[i][j] * x[j];
        }
    }
}

// ==================================================================================
// CODE2.3-XdotY.C. A C module implementing Pseudocode 2.3.                                   
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
double XdotY(int n, double x[], double y[]) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}


// ==================================================================================
// CODE3.3-CGM.C. A C module implementing Pseudocode 3.3.                                         
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
// DESCRIPTION: A module to solve Ax=b linear system with the Conjugate Gradient Method.   
//                                                                                             
// ON ENTRY                                                                                    
//    n   :: Number of equations;                                                              
//    A   :: Coefficient matrix (nÃ—n);                                                        
//    b   :: Array of length n containing the right hand side;                                 
//    x   :: Array of length n containing the initial guess;                                   
//   eps  :: Convergence tolerance;                                                            
//  maxit :: Maximum number of iterations.                                                     
//                                                                                             
// ON EXIT                                                                                     
//    x   :: Array of length n containing the approximate solution;                            
//  iter  :: Total number of iterations performed;                                             
//  error :: Euclidean (L2-) norm of displacement at exit.                                     
//                                                                                             
// USES                                                                                          
//   AX   :: A module to evaluate A*x matrix vector product;                                 
//   SQRT :: Built-in Intrinsic function returning the square root of a real value;            
//   XdotY:: A function to evaluate the dot product of two vectors.                              
//                                                                                             
// REVISION DATE :: 12/11/2024                                                                 
// ==================================================================================
void CGM(int n, double eps, double A[][4], double b[], double x[], int maxit, int* iter, double* error) {
    double c[4], r[4], d[4];
    double Enorm, rho0, rden, alpha, beta, rho;
    int p;

    AX(n, A, x, r);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i];  // [r]^(0)=[b]-[A][x]^(0)
        d[i] = r[i];         // [d]^(0)=[r]^(0)
    }
    rho0 = XdotY(n, r, r);  // rho^(0)=[r^(0)].[r^(0)]

    printf("\n\n*** Iteration history ***\n");
    for (p = 0; p <= maxit; p++) {
        Enorm = sqrt(rho0);  // E-norm=Sqrt([r]^(p-1).[r]^(p-1))
        printf("iter=%d    E-norm= %e\n", p, Enorm);
        if (Enorm < eps) break;  // Check for convergence, exit if converged

        AX(n, A, d, c);  // [c]^(p-1)=[A][d]^(p-1)
        rden = XdotY(n, d, c);  // rden=[d]^(p-1).[c]^(p-1)
        alpha = rho0 / rden;  // alpha^(p)=[r].[r]/([d].[c])
        for (int k = 0; k < n; k++) {
            x[k] = x[k] + alpha * d[k];  // [x]^(p)=[x]^(p-1)+alfa^(p)*[d]^(p-1)
            r[k] = r[k] - alpha * c[k];  // [r]^(p)=[r]^(p-1)-alfa^(p)*[d]^(p-1)
        }
        rho = XdotY(n, r, r);  // rho^(p+1)=[r^(p)].[r^(p)]
        beta = rho / rho0;  // beta^(p) = rho/rho0
        for (int k = 0; k < n; k++) {
            d[k] = r[k] + beta * d[k];  // [d]^(p)=[r]^(p)+beta*[d]^(p-1)
        }
        rho0 = rho;  // rho^(p)<==rho^(p+1)
    }

    *iter = p;  // Set total number of iterations
    *error = Enorm;  // Set recent Enorm as error
    if (*iter == maxit) {
        printf("\nFailed to converge after %d iterations\n\n", maxit);
    }
}


// ==================================================================================
// Main program to test CGM.C
// ==================================================================================
int main() {
    const int n = 4;
    double a[4][4] = {
        {3.0, -1.0, 0.0, 0.0},
        {-1.0, 3.0, -1.0, 0.0},
        {0.0, -1.0, 3.0, -1.0},
        {0.0, 0.0, -1.0, 3.0}
    };
    double b[4] = {1.5, 1.5, 2.0, 5.5};
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    double eps = 1.0e-07;
    int maxit = 999;
    int iter;
    double error;

    printf("\nLinear System\n");
    for (int i = 0; i < n; i++) {
        printf("  ");
        for (int j = 0; j < n; j++) {
            printf("%f ", a[i][j]);
        }
        printf("%f\n", b[i]);
    }

    CGM(n, eps, a, b, x, maxit, &iter, &error);

    printf("\n*** Solution ***\n");
    for (int i = 0; i < n; i++) {
        printf("  %d  %f\n", i+1, x[i]);
    }
    printf("\n Total no of iterations = %d\n Maximum Error          = %e\n", iter, error);

    return 0;
}