% ==============================================================================
%  The main program to test the module NONLINEAR_NEWTON.M
% ==============================================================================

nmax = 201;
x = zeros(nmax, 1);
y = zeros(nmax, 1);
yo = zeros(nmax, 1);
alpha = zeros(2, 1);
beta = zeros(2, 1);
gamma = zeros(2, 1);

n = 10;
if n > nmax
    fprintf('Max grid size is %d\n', nmax);
    fprintf('Increase NMAX, and then try it again.\n');
    return;
end

% Setup the ODE, Grids & BCs
eps = 1e-6;
guess = 0.4;
maxit = 99;
neq = n + 1;
xa = 0.0;
xb = 1.0;
h = (xb - xa) / n;
x = xa + (0:n)' * h;

alpha(1) = 0.0; beta(1) = 1.0; gamma(1) = 1.0;
alpha(2) = 1.0; beta(2) = 0.0; gamma(2) = 0.0;

% Prep initial guess for the solution
bc_left = abs(alpha(1));
bc_rigt = abs(alpha(2));
yo = guess * ones(neq, 1);
if bc_left == 0 && bc_rigt == 0
    yo = gamma(1) + (gamma(2) - gamma(1)) * (x - xa) / (xb - xa);
end

% Call NONLINEAR_NEWTON function
[y, ~] = NONLINEAR_NEWTON(neq, eps, x, yo, alpha, beta, gamma, maxit);

% PRINT OUT THE RESULTS
fprintf('\n%8s%8s%12s%12s\n', 'x', 'Exact', 'N.Approx', 'Abs Error');
for i = 1:neq
    error = abs(Exact(x(i)) - y(i));
    fprintf('%10.5f%12.7f%12.7f%12.3e\n', x(i), Exact(x(i)), y(i), error);
end

fprintf('\n*** D O N E ***\n\n');

function [y, error] = NONLINEAR_NEWTON(M, eps, x, yo, alpha, beta, gamma, maxit)
%  ==================================================================================
%  CODE10.4-NONLINEAR_NEWTON.m. A Matlab module implementing Pseudocode 10.4.              
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
%   DESCRIPTION: A module to find approximate solution of a two-point nonlinear differential   
%     equation using the Newton's Method. The nonlinear equations is cast in the following form
%           y'' = f(x,y,y')  on [a,b]                                                          
%     subject to                                                                               
%           alpha1 * y'(a)+ beta1 * y(a) = gamma1                                              
%           alpha2 * y'(b)+ beta2 * y(b) = gamma2                                              
%                                                                                              
%   CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1. 
%                                                                                              
%   ON ENTRY                                                                                   
%      M   :: Number of (equations) grid poinds;                                               
%      x   :: Array of length M containing the abscissas of the grid points;                   
%      yo  :: Array of length M containing the initial guess for the solution;                 
%     alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                  
%            the boundary conditions as stated above;                                          
%     eps  :: Convergence tolerance;                                                           
%     maxit:: Maximum number of iterations permitted.                                          
%                                                                                              
%   ON EXIT                                                                                    
%      y   :: Array of length M containing the approximate solution.                           
%                                                                                              
%   USES                                                                                       
%     FUNCS :: A user-defined external function module providing the coefficients of the 
%             nonlinear two-point BVP;           
%     ENORM:: A function module to calculate the Euclidean vector (L2 norm) of a vector;       
%     TRIDIAGONAL :: A module to solve a tridiagonal system of equations with Thomas algorithm.
%                                                                                              
%   REVISION DATE :: 03/10/2025                                                                
%  ==================================================================================
    a = zeros(M, 1);
    b = zeros(M, 1);
    d = zeros(M, 1);
    c = zeros(M, 1);
    y = yo;
    
    bc_left = abs(alpha(1));
    bc_rigt = abs(alpha(2));
    
    h = x(2) - x(1);
    hsqr = h^2;
    bb2h = 0.5 / h;
    hb2 = 0.5 * h;
    c1d = 2.0 * h * beta(1) / alpha(1);
    c1r = 2.0 * h * gamma(1) / alpha(1);
    c2d = 2.0 * h * beta(2) / alpha(2);
    c2r = 2.0 * h * gamma(2) / alpha(2);
    
    error = 1.0;
    p = 0;
    
    while true
        for k = 1:M
            if k == 1
                if bc_left == 0
                    d(k) = 1.0;
                    yo(k) = gamma(1) / beta(1);
                    a(k) = 0.0;
                    b(k) = 0.0;
                    c(k) = 0.0;
                else
                    yx = (gamma(1) - beta(1) * yo(k)) / alpha(1);
                    [f, fy, fp] = FUNCS(x(k), yo(k), yx);
                    a(k) = 2.0;
                    b(k) = 0.0;
                    d(k) = -2.0 + c1d - hsqr * (fy - fp * beta(1) * yo(k) / alpha(1));
                    c(k) = -c1r + (-2.0 + c1d) * yo(k) + 2.0 * yo(k+1) - hsqr * f;
                end
            elseif k == M
                if bc_rigt == 0
                    d(k) = 1.0;
                    yo(k) = gamma(2) / beta(2);
                    a(k) = 0.0;
                    b(k) = 0.0;
                    c(k) = 0.0;
                else
                    yx = (gamma(2) - beta(2) * yo(k)) / alpha(2);
                    [f, fy, fp] = FUNCS(x(k), yo(k), yx);
                    b(k) = 2.0;
                    d(k) = -2.0 - c2d - hsqr * (fy - fp * beta(2) * yo(k) / alpha(2));
                    a(k) = 0.0;
                    c(k) = c2r - (2.0 + c2d) * yo(k) + 2.0 * yo(k-1) - hsqr * f;
                end
            else
                yx = (yo(k+1) - yo(k-1)) / (2.0 * h);
                yxx = yo(k+1) - 2.0 * yo(k) + yo(k-1);
                [f, fy, fp] = FUNCS(x(k), yo(k), yx);
                b(k) = 1.0 + hb2 * fp;
                d(k) = -2.0 - hsqr * fy;
                a(k) = 1.0 - hb2 * fp;
                c(k) = yxx - hsqr * f;
            end
        end
        
        % Solve tridiagonal system of equations for correction or displacement
        c = TRIDIAGONAL(1, M, b, d, a, c);
        error = ENORM(c);
        y = yo - c;
        p = p + 1;
        fprintf('Iter=%3d   E-norm=%14.5e\n', p, error);
        yo = y;
        if error < eps || p == maxit
            break;
        end
    end
end

function [f, dfdy, dfdp] = FUNCS(x, y, yp)
% ==============================================================================
% DESCRIPTION: A user-defined function supplying the coefficients of the nonlinear 
%    two-point BVP given in the form: 
%              y'' = f(x,y,y')  on [a,b]
%
% ON ENTRY
%   x    :: Independent variable (a <=x<= b);
%   y    :: Dependent variable y=y(x).
%
% ON EXIT
%   f    :: The nonlinear two-point BVP f(x,y,y') evaluated at (x,y);
%   yp   :: First derivative of the dependent variable y';
%   dfdy :: Partial derivative of f wrt y evaluated at (x,y), df/dy;
%   dfdp :: Partial derivative of f wrt y' evaluated at (x,y), df/dy'.
%
% REVISION DATE :: 03/10/2025
% ==============================================================================
    f = y - yp^2 / y;
    dfdy = 1.0 + (yp / y)^2;
    dfdp = -2.0 * yp / y;
end

function x = TRIDIAGONAL(s1, sn, b, d, a, c)
%  ==================================================================================
%  CODE2.13-TRIDIAGONAL.m. A Matlab script module implementing Pseudocode 2.13.                   
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
%  DESCRIPTION: A subroutine to solve a tridiagonal system of linear equations                 
%    using Thomas algorithm.                                                                   
%                                                                                              
%  ON ENTRY                                                                                    
%     s1 :: Subscript of the first unknown (usually 1);                                        
%     sn :: Subscript of the last unknown (usually no. of eqs, n);                             
%      b :: Array of length n containing coefficients of below diagonal elements;              
%      d :: Array of length n containing coefficients of diagonal elements;                    
%      a :: Array of length n containing coefficients of above diagonal elements;              
%      c :: Array of length n containing coefficients of rhs.                                  
%                                                                                              
%  ON RETURN                                                                                     
%      x :: An array of length n containing the solution.                                      
%                                                                                              
%  REVISION DATE :: 03/18/2024                                                                 
%  ==================================================================================
    x = zeros(sn, 1);
    for i = s1+1:sn
        ratio = b(i) / d(i-1);
        d(i) = d(i) - ratio * a(i-1);
        c(i) = c(i) - ratio * c(i-1);
    end
    
    x(sn) = c(sn) / d(sn);
    for i = sn-1:-1:s1
        x(i) = (c(i) - a(i) * x(i+1)) / d(i);
    end
end

function norm = ENORM(x)
% ==================================================================================
% CODE3.1-ENORM.C. A C module implementing ENORM of Pseudocode 3.1.                                      
%
% NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
% First Edition. (c) By Zekeriya ALTAÇ (2024).
% ISBN: 978-1-032-75474-1 (hbk)
% ISBN: 978-1-032-75642-4 (pbk)
% ISBN: 978-1-003-47494-4 (ebk)
% 
% DOI : 10.1201/9781003474944
% C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
% 
% This free software is complimented by the author to accompany the textbook.
% E-mail: altacz@gmail.com.
% 
%  DESCRIPTION: A function module to compute Euclidean (L2) norm of a vector.                  
%                                                                                             
%  ARGUMENTS                                                                                   
%     n  :: The length of an input vector;                                                    
%     x  :: A vector (array) of length n.                                                     
%                                                                                             
%  USES                                                                                        
%   sqrt :: Built-in Intrinsic function returning the square root of a real value.            
%                                                                                             
%  REVISION DATE :: 11/09/2024                                                                  
% ==================================================================================
    norm = sqrt(sum(x.^2));
end

function y = Exact(x)
% ==============================================================================
%  DESCRIPTION: A function subprogram providing the true solution y=f(x) for 
%    testing the example problem. 
%
%  ARGUMENTS:
%     x   :: A double input, independent variable.
%
%  USES                                                                                        
%    sqrt :: Built-in intrinsic function returning the square root of a real value;
%    cosh :: Built-in Intrinsic function returning the hyperbolic cosine of a real value.  
% ==============================================================================
    aa = sqrt(2.0);
     y = sqrt(cosh(aa * (1.0 - x)) / cosh(aa));
end