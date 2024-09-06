% ==============================================================================
%  The main program to test "quadratic_eq"
% ==============================================================================
clear; clc;

% provide inputs
p = 3.; q = 2.;

% Solve quadratic equation
[xr, xi] = quadratic_eq(p, q);

% Print results
fprintf("1st Root %9.4f + %9.4f i\n", xr(1), xi(1));
fprintf("2nd Root %9.4f + %9.4f i\n", xr(2), xi(2));

function [re, im] = quadratic_eq(p, q)
% ==============================================================================
%  CODE1.3-Quadratic_Eq.m. A Mat Lab module for implementing Pseudocode 1.3.            
% 
%  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS With Pseudocodes.
%  First Edition. (c) By Zekeriya ALTAC (2024).
%  ISBN: 978-1-032-75474-1 (hbk)
%  ISBN: 978-1-032-75642-4 (pbk)
%  ISBN: 978-1-003-47494-4 (ebk)
%  
%  DOI : 10.1201/9781003474944
%  C&H/CRC PRESS, Boca Raton & London.
%  
%  This free software is complimented by the author to accompany the textbook.
%  E-mail: altacz@gmail.com.
%                                                                                    
%  DESCRIPTION: A Subroutine to find the roots of a quadratic equation of            
%    the form : x*x + p * x + q = 0.                                                 
%                                                                                    
%  INPUTS                                                                          
%    p, q :: Coefficients of the quadratic equation;                                 
%                                                                                    
%  OUTPUTS                                                                           
%    re   :: Array of length 2 containing real parts of the roots: re1, re2;         
%    im   :: Array of length 2 containing imaginary parts of the roots: im1, im2.    
%                                                                                    
%  USES                                                                              
%    sqrt :: Built-in Intrinsic function to evaluate the square root of a real value.
%                                                                                    
%  REVISION DATE :: 03/18/2024   
% ==============================================================================

d = p*p - 4*q;
if d < 0
    d = sqrt(-d);
    re = [-p/2, -p/2];
    im = [-d/2, d/2];
else
    d = sqrt(d);
    re = [(-p-d)/2, (-p+d)/2];
    im = [0, 0];
end
end