% PROGRAM Test_Back_Substitution
clear; clc;

n = 5;
a = [1.0, 2.0, 4.0, -3.0, 1.0; 
     0.0, 3.0, 4.0, -4.0, 1.0;
     0.0, 0.0, 4.0, 3.0, 1.0;
     0.0, 0.0, 0.0, 2.0, 1.0;
     0.0, 0.0, 0.0, 0.0, 1.0];

b = [10.0, 7.0, 29.0, 13.0, 5.0];

for i = 1:n
    fprintf('%10.5f %10.5f %10.5f %10.5f %10.5f %10.2f\n', a(i,1), a(i,2), a(i,3), a(i,4), a(i,5), b(i));
end
fprintf(' ********* End of Input data *********\n\n');

x = back_substitute(n, a, b);

fprintf('------ matrix A ---------------\n');
for i = 1:n
    fprintf('%10.5f %10.5f %10.5f %10.5f %10.5f\n', a(i,1), a(i,2), a(i,3), a(i,4), a(i,5));
end
fprintf('------- solution ------------------\n');
for i = 1:n
    fprintf('%d %10.5f\n', i, x(i));
end
fprintf('-----------------------------------\n');

function x = back_substitute(n, a, b)
%  ==================================================================================
%  CODE2.8-BACK_SUBSTITUTE.m. A Matlab script module implementing Pseudocode 2.8.                 
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
%  DESCRIPTION: A function module to find the solution of a upper-triangular system                 
%    using back substitution algorithm.                                                        
%                                                                                              
%  ON ENTRY                                                                                    
%      n  :: Number of unknowns;                                                               
%      a  :: Input coefficient (upper-triangular) matrix (nÃ—n);                               
%      b  :: Input array of size n containing the rhs.                                         
%                                                                                              
%  ON RETURN                                                                                     
%      x  :: Output array of size n containing the solution.                                   
%                                                                                              
%  USES
%      zeroes :: MatLab function to create array of zeroes.
%
%  REVISION DATE :: 03/18/2024                                                                 
%  ==================================================================================

x = zeros(n, 1);
x(n) = b(n) / a(n,n);
for k = n-1:-1:1
    sums = 0.0;
    for j = k+1:n
        sums = sums + a(k,j) * x(j);
    end
    x(k) = (b(k) - sums) / a(k,k);
end
end