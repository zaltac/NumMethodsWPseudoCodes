% MATLAB CODE TRANSLATED FROM FORTRAN

% ----------------------------------------------------------------
%        MATLAB APPLICATION PSEUDOCODES 2.09  
%        GAUSS ELIMINATION ON FULL SQUARE MATRIX THAT EMPLOYS
%        WITHOUT AND WITH PARTIAL PIVOTING
% ----------------------------------------------------------------        
clear; clc;

n = 7;
a = [ 2.0, 1.0, -1.0,  0.0,  3.0,  1.0, 0.0; 
     -1.0, 3.0,  1.0,  2.0,  4.0, -2.0, 1.0;
      4.0, 1.0,  5.0,  1.0, -3.0,  2.0, 2.0;
     -2.0, 1.0, -2.0, -2.0,  3.0,  3.0, 1.0;
      1.0, 1.0,  1.0,  3.0,  2.0, -3.0, 2.0;
      4.0, 1.0, -1.0,  0.0,  3.0, -5.0, 1.0;
      2.0, 7.0, -3.0, -4.0, -1.0,  2.0, 2.0];
b = [-10.0, 7.0, 28.0, -30.0, 16.0, 3.0, -12.0];
x = zeros(n, 1);

for i = 1:n
    fprintf('%8.3f ', a(i,:), b(i));
    fprintf('\n');
end

opt = 1;
% ***    SOLVE LINEAR SYSTEM OF EQS WITH NAIVE GAUSS-ELIMINATION         
x = linear_solve(n, a, b, opt);

if opt == 0
    fprintf('Naive Gauss Elimination %d\n', opt);
elseif opt == 1
    fprintf('Gauss Jordan Elimination %d\n', opt);
end

fprintf('Solution Vector\n');
for i = 1:n
    fprintf('x(%d)=%.3f\n', i, x(i));
end

function x = linear_solve(n, a, b, opt)
%  ==================================================================================
%  CODE2.9-LINEAR_SOLVE.m. A Matlab script module implementing Pseudocode 2.9.                    
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
%  DESCRIPTION: A subroutine to solve a system of linear equations using naive                 
%    Gauss Elimination (opt=0) or Gauss-Jordan Elimination (opt/=0) algorithm.                 
%                                                                                              
%  ON ENTRY                                                                                    
%      n  :: Number of unknowns;                                                               
%      A  :: Input coefficient matrix of size nÃ—n;                                            
%      b  :: An input array of length n containing the rhs;                                    
%     opt :: Option key (=0, Naive Gauss-Elim.; /=0, Gauss-Jordan Elimn).                      
%                                                                                              
%  ON RETURN                                                                                     
%      x  :: The output array of length n containing the solution.                             
%                                                                                              
%  USES                                                                                                                                                                           
%    abs  :: Built-in Intrinsic function returning the absolute value of a real value;         
%    back_substitute :: A subrotine to solve an upper-triangular system.                       
%                                                                                              
%  REVISION DATE :: 03/18/2024                                                                 
%  ==================================================================================
eps = 1e-12;

%  *** FORWARD ELIMINATION STEPS
for j = 1:n
    ajj = a(j,j);
    if abs(ajj) < eps
        fprintf('Pivot is zero at j=%d\n', j);
        fprintf('Execution is halted!\n');
        return;
    else
        val = 1/ajj;
        b(j) = b(j)*val;
        a(j,:) = a(j,:)*val;
        for i = j+1:n
            s = a(i,j);
            a(i,j) = 0;
            a(i,j+1:n) = a(i,j+1:n) - s*a(j,j+1:n);
            b(i) = b(i) - s*b(j);
        end
    end
end

% *** BACK SUBSTITUTION STEPS
if opt == 0
    x = back_substitute(n, a, b);
else
    for j = n:-1:2
        for i = j-1:-1:1
            b(i) = b(i) - a(i,j)*b(j);
            a(i,j) = 0;
        end
        x(j) = b(j);
    end
    x(1) = b(1);
end
end

function x = back_substitute(n, a, b)
% --------------------------------------------------------------------
% DESCRIPTION: SOLVES UPPER TRIANGULAR SYSTEM USING FORWARD SUBSTITUTION
%   (applies Pseudocodes 2.8)
%
% On ENTRY
%    n :: Number of unknowns; 
%    A :: Input coefficient matrix of size n×n;
%    B :: An array of size n containing the rhs.
%     
% On RETURN
%    X  :: An array of size n containing the solution.
%
%  USES                                                                                             
%    zeros :: Built-in function creating an mxn matrix of zeros. 
% --------------------------------------------------------------------
x = zeros(n, 1);
x(n) = b(n)/a(n,n);
for k = n-1:-1:1
    sums = 0;
    for j = k+1:n
        sums = sums + a(k,j)*x(j);
    end
    x(k) = (b(k) - sums)/a(k,k);
end
end