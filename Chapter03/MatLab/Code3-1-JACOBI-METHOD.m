% ==============================================================================
%  The main program to test JACOBI.M
% ==============================================================================
% Import necessary libraries
import matlab.io.*

% Define constants
n = 10;

% Initialize matrices and vectors
a = [2 -1  0  0  0  0  0  0  0  0;
    -1  2 -1  0  0  0  0  0  0  0;
     0 -1  2 -1  0  0  0  0  0  0;
     0  0 -1  2 -1  0  0  0  0  0;
     0  0  0 -1  2 -1  0  0  0  0;
     0  0  0  0 -1  2 -1  0  0  0;
     0  0  0  0  0 -1  2 -1  0  0;
     0  0  0  0  0  0 -1  2 -1  0;
     0  0  0  0  0  0  0 -1  2 -1;
     0  0  0  0  0  0  0  0 -1  2];

b = [0 0 0 0 0 0 0 0 0 2.2]';
xo = zeros(n, 1);

% Print input data
fprintf('    Coefficient Matrix    RHS    Initial Guess\n');
for i = 1:n
    fprintf('%s %12.2f %12.2f\n', sprintf('%6.1f', a(i,:)), b(i), xo(i));
end

% Set parameters
eps = 1e-7;
maxit = 999;
Iprint = 1;
iter = Iprint;

% Call Jacobi method
[x, iter, Errmax] = jacobi(n, eps, a, b, xo, maxit);

% Print results
fprintf('\n----- Solution --------\n');
for i = 1:n
    fprintf('%3d  %11.6f\n', i, x(i));
end
fprintf('-----------------------\n');
fprintf('Total no of iterations = %4d\n', iter);
fprintf('Maximum Error = %10.3e\n\n', Errmax);



%  ==================================================================================
%  CODE3.1-JACOBI.m. A Matlab script module implementing Pseudocode 3.1.                          
% 
%  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
%  First Edition. (c) By Zekeriya ALTAÇ (2024).
%  ISBN: 978-1-032-75474-1 (hbk)
%  ISBN: 978-1-032-75642-4 (pbk)
%  ISBN: 978-1-003-47494-4 (ebk)
%  
%  DOI : 10.1201/9781003474944
%  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
%  
%  This free software is complimented by the author to accompany the textbook.
%  E-mail: altacz@gmail.com.
%  
%  DESCRIPTION: A function to iteratively solve Ax=b using the Jacobi method.                
%                                                                                              
%  ON ENTRY                                                                                    
%     n   :: Number of equations (size of A);                                                  
%     A   :: Input coefficient matrix (nÃ—n);                                                  
%     b   :: Array of length n containing the right-hand;                                      
%     x   :: Array of length n containing the estimate at (p+1)'th step;                       
%     xo  :: Array of length n containing the initial guess, or iterates at estimate at p'th st
%    eps  :: Convergence tolerance;                                                            
%   maxit :: Maximum permitted number of iterations.                                           
%                                                                                              
%  ON EXIT                                                                                     
%     x   :: Array of length n containing the estimated solution;                              
%   iter  :: Total number of iterations performed;                                             
%   error :: L2 norm of the displacement vector.                                               
%                                                                                              
%  USES                                                                                        
%    JACOBI_DRV :: Accompanying function performing one step Jacobi iteration.               
%                                                                                              
%  REVISION DATE :: 11/09/2024                                                                 
%  ==================================================================================
function [x, iter, error] = jacobi(n, eps, A, b, xo, maxit)
    x = xo;
    del0 = 1.0;
    
    for p = 1:maxit
        [x, del] = jacobi_drv(n, A, b, xo);
        
        % Print iteration progress
        fprintf('iter=%4d    Error= %11.5e    Ratio= %11.5e\n', p, del, del/del0);
        
        xo = x;
        del0 = del;
        
        if del < eps
            break;
        end
    end
    
    error = del;
    iter = p;
    
    if p == maxit
        fprintf('\nJacobi method failed to converge after %d iterations\n', maxit);
        fprintf('within the specified EPS tolerance.\n\n');
    end
end


%  ==================================================================================
%  CODE3.1-JACOBI.m. A Matlab script module implementing JACOBI_DRV in Pseudocode 3.1.                          
% 
%  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
%  First Edition. (c) By Zekeriya ALTAÇ (2024).
%  ISBN: 978-1-032-75474-1 (hbk)
%  ISBN: 978-1-032-75642-4 (pbk)
%  ISBN: 978-1-003-47494-4 (ebk)
%  
%  DOI : 10.1201/9781003474944
%  C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
%  
%  This free software is complimented by the author to accompany the textbook.
%  E-mail: altacz@gmail.com.
%  
%  DESCRIPTION: A function module to perform one step Jacobi iteration and compute               
%     the Euclidean norm of the displacement vector.                                           
%                                                                                              
%  ON ENTRY                                                                                    
%     n   :: Number of equations (size of A);                                                  
%     A   :: Input coefficient matrix (nÃ—n);                                                  
%     b   :: Array of length n containing the right-hand;                                      
%     x   :: Array of length n containing the estimate at (p+1)'th step;                       
%     xo  :: Array of length n containing the estimate at p'th step.                           
%                                                                                              
%  ON EXIT                                                                                     
%     x   :: Array of length n containing the estimated solution;                              
%     del :: Maximum absolute error achieved.                                                  
%                                                                                              
%  USES                                                                                        
%   NORM:: Built-in function calculating the Euclidean (L2) norm of a vector.      
%                                                                                              
%  REVISION DATE :: 11/09/2024                                                                  
%  ==================================================================================
function [x, del] = jacobi_drv(n, A, b, xo)
    x = zeros(n, 1);
    d = zeros(n, 1);
    
    for i = 1:n
        sums = 0;
        for j = 1:n
            if i ~= j
                sums = sums + A(i,j) * xo(j);
            end
        end
        x(i) = (b(i) - sums) / A(i,i);
        d(i) = x(i) - xo(i);
    end
    
    del = norm(d);
end