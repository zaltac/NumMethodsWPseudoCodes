% ==============================================================================
%  The main program to test LBVP_SOLVE.M
% ==============================================================================
nmax = 501;

% Initialize variables
x = zeros(1, nmax);
y = zeros(1, nmax);
alpha = zeros(1, 2);
beta = zeros(1, 2);
gamma = zeros(1, 2);

% Input number of intervals
n = 10;
if n > nmax
    fprintf('Max grid size is %d\n', nmax);
    fprintf('Increase NMAX, and then try it again.\n');
    return;
end

% Setup the ODE, Grids & BCs
neq = n + 1;
a = 1.0;
b = 2.0;
h = (b - a) / n;
x(1) = a;
for i = 1:n
    x(i+1) = x(i) + h;
end

alpha(1) = 1.0; beta(1) = 2.0; gamma(1) = 4.0;
alpha(2) = 1.0; beta(2) =-2.0; gamma(2) = 2.0;

if (abs(alpha(1)) > 0.0 && alpha(1) < 1.0) || (abs(alpha(2)) > 0.0 && alpha(2) < 1.0)
   fprintf(' ******* WARNING alfa should 0 or >=1 *********\n');
end

% Solve the BVP
y = LBVP_SOLVE(neq, x, alpha, beta, gamma);

% Print results
fprintf('%8s%7s%12s%16s\n', 'x', 'Exact', 'N.Approx', 'Abs Error');
for i = 1:neq
    error = abs(y(i) - exact(x(i)));
    fprintf('%10.4f%12.7f%12.7f%16.5e\n', x(i), exact(x(i)), y(i), error);
end


function y = LBVP_SOLVE(M, x, alpha, beta, gamma)
%  ==================================================================================
%  CODE10.1-LBVP_SOLVE.m. A Matlab script module implementing Pseudocode 10.1.                    
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
%  DESCRIPTION: A module to find approximate solution of the following linear                  
%    differential equation using the Finite Difference Method:                                 
%          p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]                                   
%    subject to                                                                                
%          alpha1 * y'(a)+ beta1 * y(a) = gamma1                                               
%          alpha2 * y'(b)+ beta2 * y(b) = gamma2                                               
%                                                                                              
%  CAUTION!!! In case of alpha<>0 make sure that the BCs are normalized so that alphas are 1.  
%                                                                                              
%  ON ENTRY                                                                                    
%    neq  :: Number of (equations) grid poinds;                                                
%     x   :: Array of length neq containing the abscissa of the grid points;                   
%    alpha, beta, gamma :: Arrays of length 2 containing the coefficients of                   
%           the boundary conditions as stated above;                                           
%                                                                                              
%  ON RETURN                                                                                     
%     y   :: Array of length neq containing the approximate solution.                          
%                                                                                              
%  USES                                                                                        
%   COEFFS  :: A module containing the coefficients and rhs of the linear ODE, i.e., p(x), q(x)
%   TRIDIAGONAL:: A module solving a tridiagonal system of equations using the Thomas algorithm
%                                                                                              
%  REVISION DATE :: 03/08/2025                                                                 
%  ==================================================================================
    dx = diff(x);
    C = zeros(1, M);
    y = zeros(1, M);
    B = zeros(1, M);
    A = zeros(1, M);
    D = zeros(1, M);

    BC_Left = floor(alpha(1));
    BC_Rigt = floor(alpha(2));

    for k = 1:M
        xk = x(k);
        hp = dx(min(k, M-1));
        xp = xk + 0.5 * hp;
        [pp, qp, rp, gp] = COEFFS(xp);

        if k == 1  % First half-control volume
            if BC_Left == 0
                D(k) = 1.0;
                A(k) = 0.0;
                B(k) = 0.0;
                C(k) = gamma(1) / beta(1);
            else
                A(k) = pp / hp + 0.5 * qp;
                B(k) = 0.0;
                D(k) = -A(k) + 0.5 * rp * hp + pp * beta(1) / alpha(1);
                C(k) = 0.5 * gp * hp - pp * gamma(1) / alpha(1);
            end
        elseif k == M  % Last half-control volume
            if BC_Rigt == 0
                D(k) = 1.0;
                A(k) = 0.0;
                B(k) = 0.0;
                C(k) = gamma(2) / beta(2);
            else
                B(k) = pm / hm - 0.5 * qm;
                A(k) = 0.0;
                D(k) = -B(k) - pm * beta(2) / alpha(2) + 0.5 * rm * hm;
                C(k) = 0.5 * gm * hm - pm * gamma(2) / alpha(2);
            end
        else  % Interior control volumes
            B(k) = pm / hm - 0.5 * qm;
            A(k) = pp / hp + 0.5 * qp;
            D(k) = -A(k) - B(k) + 0.5 * (rp * hp + rm * hm);
            C(k) = 0.5 * (gp * hp + gm * hm);
        end

        hm = hp;
        pm = pp; qm = qp;
        rm = rp; gm = gp;
    end

    y = TRIDIAGONAL(1, M, B, D, A, C);
end

function [p, q, r, g] = COEFFS(x)
% ==============================================================================
% DESCRIPTION: A user-defined module suppling the coefficients p(x), q(x), r(x) and
%    rhs g(x) of the linear ODE given in the following form: 
%         p(x) * y''+ q(x) * y' + r(x) * y = g(x)  on [a,b]
%
% ON ENTRY
%    x   :: Independent variable (a<= x<= b);
%
% ON RETURN
%   p, q, r, abd g :: Coefficients & rhs evaluated at x. 
%
% REVISION DATE :: 03/08/2025
% ==============================================================================
    p = x*x;
    q = -5.0*x;
    r = 8.0;
    g = 0.0;
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
    x = zeros(1, sn);
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

function y = exact(x)
% ==============================================================================
% DESCRIPTION: A function subprogram providing the true solution y=f(x) for testing the module. 
%
% ARGUMENTS:
%      x   :: A real input, independent variable.
% ==============================================================================
    y = x^2*(x^2 - 0.5);
end

 
