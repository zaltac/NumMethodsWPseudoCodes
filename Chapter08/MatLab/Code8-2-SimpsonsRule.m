% ==============================================================================
%  The main program to test function module Simpsons_Rule_DF.M
% ==============================================================================
clear; clc;

n = input("Enter Number of Panels: ");
a = 0.0;
b = 1.0;
h = (b-a)/n;
x = a;
f = zeros(1,n+1);
for i = 1:n+1   % *** Construct a discrete function using f(x)=x^4
    f(i) =x^4;
    x = x + h;
end

[intg] = Simpsons_Rule_DF(n,h,f);

fprintf("=== Standard Simpsons Rule ===\n");
fprintf("%d %f\n\n", n, intg);



%  ==================================================================================
%  CODE8.2-Simpsons_Rule_DF.m. A Matlab script module implementing Pseudocode 8.2.                
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
%  DESCRIPTION: A function module to estimate the integral of a discrete function f on [a,b]                         
%   using the Simpson's 1/3 rule.                                                              
%                                                                                              
%  ON ENTRY                                                                                    
%     n   :: Number of panels (must be even!..);                                               
%     h   :: Interval size (uniform spacing, x_{i+1}-x_i);                                     
%     f   :: Array of length (n+1) containing the ordinates, f_1, f_2, ..., f_{n+1}.               
%                                                                                              
%  ON EXIT                                                                                     
%   intg  :: Numerical estimate for the integral.                                              
%                                                                                              
%  USES                                                                                        
%    MOD  :: Modulo function, returning the remainder after number is divided by divisor.      
%                                                                                              
%  REVISION DATE :: 03/04/2025                                                                 
%  ==================================================================================
function [intg] = Simpsons_Rule_DF(n,h,f) 
    m = mod(n, 2);
    if m ~= 0
        error('Number of panels is not EVEN');
    end

    odd = 0.0;
    for i = 2:2:n
        odd = odd + f(i);
    end

    even = 0.0;
    for i = 3:2:n-1
        even = even + f(i);
    end

    intg = f(1) + f(n + 1) + 4.0 * odd + 2.0 * even;
    intg = intg * h / 3.0;
end