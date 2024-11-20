% PROGRAM test_Newton_Raphson
% ==============================================================================
%  The main program to test function NEWTON_RAPHSON
% ==============================================================================
clear; clc;

maxit = 99;
eps = 0.50e-4;
root = 2.20;

[root, iter] = newton_raphson(root, maxit, eps);

fprintf('-----------------------------\n');
fprintf('Root is %f, converged after %d iterations\n', root, iter);
fprintf('-----------------------------\n');



%  ==================================================================================
%  CODE4.3-NEWTON_RAPHSON.M. A Matlab script module implementing Pseudocode 4.3.                  
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
%  DESCRIPTION: A Matlab module to compute a root of a nonlinear equation using the                   
%    Newton-Raphson method.                                                                    
%                                                                                              
%  ON ENTRY                                                                                    
%    root  :: Initial guess for the root;                                                      
%    maxit :: Maximum number of iterations permitted;                                          
%    eps   :: Convergence tolerance.                                                           
%                                                                                              
%  ON RETURN                                                                                     
%    iter  :: Number of iterations realized;                                                   
%    root  :: Computed approximation for the root.                                             
%                                                                                              
%  USES                                                                                        
%    ABS   :: Built-in Intrinsic function returning the absolute value of a real value;        
%                                                                                              
%  ALSO REQUIRED                                                                               
%    FUNC  :: User-defined external function providing the nonlinear equation, f(x).           
%    FUNCP :: User-defined external function providing the first derivative                    
%             of the nonlinear equation, f'(x).                                                
%                                                                                              
%  REVISION DATE :: 11/20/2024                                                                 
%  ==================================================================================
function [root, iter] = newton_raphson(root, maxit, eps)

fprintf('   p     x^(p)       f(x^(p))     f''(x^(p))     Aerr     Rate\n');

del0 = 1.0;
x0 = root;   % Initialize root 
p = 0;       % Initialize iteration counter

% BEGIN REPEAT-UNTIL LOOP
while true
    fn = func(x0);
    fpn = funcp(x0);
    del = -fn/fpn;
    aerr = abs(del);      % absolute error
    rate = aerr/del0^2;   % estimate the convergence rate
    fprintf('%4d %12.7f %12.7f %12.7f %10.7f %10.7f\n', p, x0, fn, fpn, aerr, rate);
    xn = x0 + del;        % update the estimate
    x0 = xn;
    del0 = abs(del);
    p = p + 1;
    if (abs(fn) < eps && aerr < eps) || (p == maxit)
        break; % End the iteration loop
    end
end

root = xn;
iter = p;

if p == maxit
    fprintf('** Max iteration number reached=\n');
    fprintf('** Estimated root has NOT converged, del, f(x)= %f, %f\n', abs(del), func(xn));
end

end

function f = func(x)
% ==========================================================================          
% User-defined function providing f(x), which should be cast as func(x)=0.
% ========================================================================== 
    f = 4.0 + x^2 * (8.0 - x^2);
end

function fp = funcp(x)
% ==========================================================================          
% User-defined function providing f'(x)=funcp(x).
% ========================================================================== 
    fp = 4.0 * x * (4.0 - x^2);
end