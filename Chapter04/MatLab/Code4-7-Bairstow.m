% Test_Bairstow
% ==============================================================================
%  Main program to test BAIRSTOW function
% ==============================================================================

% Initialize variables
A = zeros(1, 20);
XRE = zeros(1, 19);
XIM = zeros(1, 19);
n = 5;
A(1:6) = [1.0, -5.0, -15.0, 85.0, -26.0, -120.0];

% Output control key
iprnt = 2;
% iprnt = 0  does not print iteration details, 
%       = 1  prints a short iteration history 
%       = 2  prints all iteration history  

maxit = 99;
p0 = 0.0;
q0 = 0.0;
eps = 0.5e-4;

%  ==================================================================================
%  CODE4.7-BAIRSTOW.M. A Matlab script module implementing Pseudocode 4.7.                        
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
%  DESCRIPTION: A function to find all real and/or imaginary roots of a polynomial           
%    of the n'th degree using the BAIRSTOW's method.                                         
%                                                                                              
%  ON ENTRY                                                                                    
%    n    :: Degree of the polynomial;                                                         
%   p0,q0 :: Initial guesses for a quadratic equation; i.e., for p and q;                      
%     a   :: Array of length (n+1) containing the coefficients of polynomial defined as        
%                 a0 x^n + a1 x^(n-1) + ... + an = 0                                           
%    eps  :: Convergence tolerance;                                                            
%   maxit :: Maximum number of iterations permitted;                                           
%   iprnt :: printing key, =0 do not print intermediate results, <> 0 print intermediates.     
%                                                                                              
%  ON RETURN                                                                                     
%    xre  :: Array of length n containing real parts of the roots;                             
%    xim  :: Array of length n containing imaginary parts of the roots.                        
%                                                                                              
%  OTHER VARIABLES                                                                             
%     b   :: Array of length [n] containing coefficients of quotient polynomial (0<=k<=n-2);   
%     c   :: Array of length [n] containing coefficients of partial derivatives.               
%                                                                                              
%  USES                                                                                        
%    abs  :: Built-in Intrinsic function returning the absolute value of a real value;         
%    QUADRATIC :: Subroutine that solves a quadratic equation of the form x2 + p x + q = 0. (see CODE1-3)   
%                                                                                              
%  REVISION DATE :: 04/29/2024                                                                 
%  ==================================================================================
[XRE, XIM] = BAIRSTOW(n, p0, q0, A, eps, maxit, iprnt);

% Print results
disp('    ======== All the Roots are =========');
for i = 1:n
    fprintf('    Root(%2d) = %8.5f + ( %8.5f ) i\n', i, XRE(i), XIM(i));
end
disp('    ====================================');

function [xre, xim] = BAIRSTOW(n, p0, q0, a, eps, maxit, iprnt)
    b = zeros(1, n+1);
    c = zeros(1, n+1);
    xre = zeros(1, n);
    xim = zeros(1, n);
    xr = zeros(1, 2);
    xi = zeros(1, 2);

    % Normalize a's by a(1)
    a = a / a(1);

    m = n; % Save n for later use
    kount = 0;
    while n > 1
        p = p0;
        q = q0; % Initialize
        k = 0;
        delM = 1.0;
        while delM > eps && k <= maxit % Inner loop
            k = k + 1;
            b(1) = 1.0;
            c(1) = 1.0;
            b(2) = a(2) - p;
            c(2) = b(2) - p;
            for i = 3:n+1
                b(i) = a(i) - p*b(i-1) - q*b(i-2);
                c(i) = b(i) - p*c(i-1) - q*c(i-2);
            end
            cbar = c(n) - b(n);
            del = c(n-1)^2 - cbar*c(n-2);
            del1 = b(n)*c(n-1) - b(n+1)*c(n-2);
            del2 = b(n+1)*c(n-1) - b(n)*cbar;
            delp = del1/del;
            delq = del2/del;
            p = p + delp;
            q = q + delq; % Find new estimates
            delM = abs(delp) + abs(delq); % Calculate L1 norm
            if iprnt == 1
                fprintf('Iter= %2d  delM= %10.4e  p= %12.5e  q= %12.5e\n', k, delM, p, q);
            elseif iprnt == 2
                fprintf('\n    iter=%4d\n    ---------\n', k);
                fprintf('     dp =%14.6e      dq =%14.6e      delM=%14.6e\n', delp, delq, delM);
                fprintf('      p =%14.6e       q =%14.6e\n\n', p, q);
                fprintf('     k          a(k)           b(k)           c(k)\n');
                fprintf('    %s\n', repmat('-', 1, 47));
                for i = 1:n+1
                    fprintf('    %2d   %12.6f   %12.6f   %12.6f\n', i-1, a(i), b(i), c(i));
                end
                fprintf('    %s\n', repmat('-', 1, 47));
            end
        end
        if k-1 == maxit
            fprintf('Quadratic factor did not converge after %d iterations\n', k-1);
            fprintf('Recent values of p, q, delM are %f, %f, %f\n', p, q, delM);
            fprintf('Corresponding roots may be questionable ...\n');
        end
        [xr, xi] = QUADRATIC_EQ(p, q);
        kount = kount + 1;
        xre(kount) = xr(1);
        xim(kount) = xi(1);
        kount = kount + 1;
        xre(kount) = xr(2);
        xim(kount) = xi(2);

        fprintf('        ======== FOUND A QUADRATIC FACTOR  ======== \n');
        fprintf('          x*x + (%12.6f)*x + (%12.6f) \n', p, q);
        fprintf('        =========================================== \n\n\n');

        n = n - 2;
        a(1:n+1) = b(1:n+1);
        if n == 1
            kount = kount + 1;
            xre(kount) = -a(2);
            xim(kount) = 0.0;
            fprintf('        ======== FOUND A LINEAR FACTOR  ======== \n');
            fprintf('               x  + (%12.6f) \n', a(2));
            fprintf('        ======================================== \n\n\n');
        end
    end
    n = m;
end

 