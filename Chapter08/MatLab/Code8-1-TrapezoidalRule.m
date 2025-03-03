% ==============================================================================
%  The main program to test Trapezoidal_Rule_RF.M
& ==============================================================================
clear; clc;
a = 0.0;
b = 1.0;
fprintf("Enter Number of Panels ");
n = input("");

[intg, intgc] = Trapezoidal_Rule_RF(n, a, b);

fprintf("=== Standard Trapezodial Rule ===\n");
fprintf("%d %f\n", n, intg);
fprintf(" \n");
fprintf("=== Trapezodial Rule with End-Correction ===\n");
fprintf("%d %f\n", n, intgc);
fprintf(" \n");


%  ==================================================================================
%  CODE8.1-Trapezoidal_Rule_RF.m. A Matlab function implementing Pseudocode 8.1.             
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
%  DESCRIPTION: A function to estimate the integral of y=f(x) on [a,b]                       
%   using the Trapezoidal rule with/without end correction.                                    
%                                                                                              
%  ON ENTRY                                                                                    
%     n   :: Number of panels (i.e., n+1 integration points);                                  
%   [a,b] :: Integration interval.                                                             
%                                                                                              
%  ON RETURN                                                                                     
%   intg   :: Integral estimate using the ordinary Trapezoidal rule;                           
%   intgc  :: Integral estimate using the Trapezoidal rule with the end-point correction.      
%                                                                                              
%  USES                                                                               
%     FX   :: User-defined external function providing the function, f(x);                     
%     FU   :: User-defined external function providing the first derivative, f'(x).            
%                                                                                              
%  REVISION DATE :: 03/03/2025                                                                 
%  ==================================================================================
function [intg, intgc] = Trapezoidal_Rule_RF(n, a, b)
h = (b - a) / n;
intg = 0.5 * (FX(a) + FX(b));
xi = a;
for i = 1:n-1
   xi = xi + h;
   intg = intg + FX(xi);
end
intg = h * intg;

corr = -h^2 * (FU(b) - FU(a)) / 12;
intgc = intg + corr;
end

function y = FX(x)
% ==============================================================================
% DESCRIPTION: User-defined function providing y=f(x) to be integrated. 
%
% ARGUMENTS:
%      x   :: a real input value.
% ==============================================================================
y = x^4;
end

function y = FU(x)
% ==============================================================================
% DESCRIPTION: User-defined function providing first derivative, f'(x), explicitly. 
%
% ARGUMENTS:
%      x   :: a real input value.
% ==============================================================================
y = 4 * x^3;
end



