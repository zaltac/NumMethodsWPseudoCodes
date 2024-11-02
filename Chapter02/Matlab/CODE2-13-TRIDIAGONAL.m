% PROGRAM Test_Tridiagonal
% ==============================================================================
%  The main program to test SUBROUTINE TRIDIAGONAL
% ==============================================================================              
clear; clc;
n = 9;
b = [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
d = [-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0];
a = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
c = [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -14.0];
x = zeros(n, 1);

fprintf('\n ********** Input Tridiagonal Matrix & RHS ***********\n');
for i = 1:n
   fprintf('%10.5f %10.5f %10.5f %10.5f\n', b(i), d(i), a(i), c(i));
end

s1 = 1;
sn = n;
x = TRIDIAGONAL(s1, sn, b, d, a, c);

% ***     PRINT OUT THE RESULTS
%         
fprintf('\n*** Solution ***\n\n');
for i = 1:n
   fprintf('   x(%d)=%6.3f\n', i, x(i));
end
fprintf('\n');


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
function x = TRIDIAGONAL(s1, sn, b, d, a, c)
    for i = s1+1:sn
        ratio = b(i)/d(i-1);
        d(i) = d(i) - ratio * a(i-1);
        c(i) = c(i) - ratio * c(i-1);
    end

    x(sn) = c(sn) / d(sn);
    for i = sn-1:-1:s1
        x(i) = (c(i) - a(i) * x(i+1)) / d(i);
    end
end