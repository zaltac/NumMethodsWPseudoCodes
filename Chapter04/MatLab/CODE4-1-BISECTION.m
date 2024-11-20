% Test_Bisection
% ==============================================================================
%  Main program to test function BISECTION.M
% ==============================================================================

maxit = 99;
eps = 0.50e-4;
a = 0.0;
b = 4.0;

fa = func(a);
fb = func(b);
if fa*fb > 0
    disp('No root in interval (a,b). Change the interval.');
    return;
end

[root, halves] = bisection(a, b, maxit, eps);

fprintf('%-44s\n', repmat('=', 1, 44));
fprintf('! Root is %14.8e after %3d bisections !\n', root, halves);
fprintf('%-44s\n', repmat('=', 1, 44));

root = (a*fb - b*fa) / (fb - fa);
fprintf('\n*** Root (with linear interpolation)   = %14.9e\n', root);
fprintf('*** Closeness to the root, Abs[f(root)]= %10.4e\n\n', abs(func(root)));


%  ==================================================================================
%  CODE4.1-BISECTION.M. A Matlab function module implementing Pseudocode 4.1.                       
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
%  DESCRIPTION: A function to find a root of a nonlinear equation in [a,b]                       
%    using the Bisection method.                                                               
%                                                                                              
%  ON ENTRY                                                                                    
%   [a,b] :: Initial search interval (it must bracket one root);                               
%   maxit :: Maximum number of iterations permitted;                                           
%   eps   :: Convergence tolerance.                                                            
%                                                                                              
%  ON RETURN                                                                                     
%   halves:: Number of halves realized;                                                        
%   root  :: Computed approximation for the root.                                              
%                                                                                              
%  USES                                                                                        
%   ABS   :: Built-in Intrinsic function returning the absolute value of a real value;         
%                                                                                              
%  ALSO REQUIRED                                                                               
%   FUNC  :: User-defined external function providing the nonlinear equation.                  
%                                                                                              
%  REVISION DATE :: 11/20/2024                                                                 
%  ==================================================================================
function [root, halves] = bisection(a, b, maxit, eps)
    p = 0;
    interval = b - a;
    fa = func(a);
    fb = func(b);

    fprintf('%3s %12s %12s %12s %12s %12s %11s %11s\n', 'p', 'a', 'b', 'f(a)', 'f(b)', 'xm', 'f(xm)', 'interval');
    fprintf('%s\n', repmat('-', 1, 97));

    while true
        p = p + 1;
        xm = 0.5 * (a + b);
        fm = func(xm);
        fprintf('%3d %12.7e %12.7e %12.4e %12.4e %12.7e %11.4e %11.4e\n', p, a, b, fa, fb, xm, fm, interval);
        
        if fa * fm > 0
            a = xm;
            fa = fm;
        else
            b = xm;
            fb = fm;
        end
        
        interval = 0.5 * interval;
        
        if ((abs(fm) < eps && interval < eps) || p == maxit)
            break;
        end
    end

    root = xm;
    halves = p;

    if p == maxit
        fprintf('%s\n', repmat('!', 1, 37));
        fprintf('Max iteration number reached=%3d%s\n', p, repmat('!', 1, 37));
    end
end

function y = func(x)
    % =======================================================================          
    % User-defined Function providing f(x) which should be cast as func(x)=0.
    % =======================================================================  
    y = x^2 + 0.025*x - 4.0;
end
