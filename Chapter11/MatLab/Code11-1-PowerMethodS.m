% ==============================================================================
%  The main program to test POWER_METHOD_S.m
% ==============================================================================
clear; clc;

n = 4;
A = [8 9 10 9; -13 -12 -12 -11; -18 -9 -20 -9; 11 1 10 0];
x = zeros(n, 1);
x(1) = 1;
lambda = 1.0;
maxit = 299;
eps = 1e-3;

% *** Apply Power Method  
[lambda, x, error] = POWER_METHOD_S(A, x, lambda, eps, maxit);

% *** Print out the results
if error < eps
    fprintf('\n%s\n', repmat('=', 1, 20));
    fprintf('Largest lambda value = %12.7f\n', lambda);
    fprintf('Eigenvecto    = %12.4f %12.4f %12.4f %12.4f\n', x);
    fprintf('%s\n', repmat('=', 1, 20));
else
    fprintf('Max number of iterations is reached.\n');
end


function [lambda, x, error] = POWER_METHOD_S(A, x, lambda, eps, maxit)
%  ==================================================================================
%  CODE11.1-POWER_METHOD_S.m. A Matlab script module implementing Pseudocode 11.1.                
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
%   DESCRIPTION: A MatLab module to find the dominant eigenvalue (largest in absolute value)          
%      using the Power Method with scaling technique.                                          
%                                                                                              
%   ON ENTRY                                                                                   
%      n     :: Size of the matrix;                                                            
%      A     :: A real square matrix (nxn);                                                    
%      x     :: Array of length n containing the initial guess for the eigenvector;            
%     lambda :: An initial guess for the dominant eigenvalue;                                  
%      eps   :: Convergence tolerance;                                                         
%     maxit  :: Maximum number of iterations permitted.                                        
%                                                                                              
%   ON EXIT                                                                                    
%     lambda :: Estimated dominant (largest in absolute value) eigenvalue;                     
%      x     :: Array of length n containing the estimated eigenvector;                        
%     error  :: Error, max of both L2-norm of displacement vector and relative error           
%               for eigenvalue.                                                               
%                                                                                              
%   USES                                                                                       
%     abs  :: Built-in intrinsic function returning the absolute value of a real value.  
%     norm :: Built-in intrinsic function returning the E-norm of a real array.             
%     max  :: Built-in intrinsic function returning the maximum of the input values.         
%     MAX_SIZE:: A function module providing largest (in absolute) value of a vector.          
%                                                                                              
%   REVISION DATE :: 03/14/2025                                                                
%  ==================================================================================
    n = size(A, 1);
    xn = x;
    lambdao = lambda;
    error = MAX_SIZE(n,x);
    p = 0;
    err2 = 1;    
    while error > eps && p < maxit
        p = p + 1;
        xn = A * x;
        lambda = MAX_SIZE(n,xn);
        xn = xn / lambda;
        err1 = norm(abs(xn - x));
        err2 = abs(1 - lambdao / lambda);
        error = max(err1, err2);
        fprintf('%3d %12.7f %14.6e\n', p, lambda, err2);
        lambdao = lambda;
        x = xn;
    end
end


function xmax = MAX_SIZE(n,x)
%  ==================================================================================
%  CODE11.1-POWER_METHOD_S.m of a module in Pseudocode 11.1.                
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
%   DESCRIPTION: A function module to find the largest element (in absolute value) of an array.
%                                                                                              
%   ARGUMENTS                                                                                  
%      n   :: Length of the array;                                                          
%      x   :: Array of length n.                                                            
%                                                                                              
%   USES                                                                                       
%     abs  :: Built-in intrinsic function returning the absolute value of a real value;       
%                                                                                              
%   REVISION DATE :: 03/10/2025                                                                
%  ==================================================================================
    xmax = x(1);
    for i = 2:n
        if abs(x(i)) > abs(xmax) 
            xmax = x(i);
        end
    end
end
