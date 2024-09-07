% PROGRAM Test_EXPE
% ==============================================================================
%  The main program to test FUNCTION EXPE
% ==============================================================================
clear; clc;

x = input("Enter x: ");
n = input("Enter n: ");

fprintf('e^%0.3f = %10.7f\n', x, expe(x,n));

function f = expe(x,n)
%  ==================================================================================
%  CODE1.7-ExpE.m. A Matlab script module implementing Pseudocode 1.7.                            
% 
%  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
%  First Edition. (c) By Zekeriya ALTAÇ (2024).
%  ISBN: 978-1-032-75474-1 (hbk)
%  ISBN: 978-1-032-75642-4 (pbk)
%  ISBN: 978-1-003-47494-4 (ebk)
%  
%  DOI : 10.1201/9781003474944. C&H/CRC PRESS, Boca Raton, FL, USA & London, UK.
%  
%  This free software is complimented by the author to accompany the textbook.
%  E-mail: altacz@gmail.com.
%                                                                                              
%  DESCRIPTION: A function to compute e^x using the MacLaurin series with specified            
%     number of terms.                                                                         
%                                                                                              
%  ARGUMENTS                                                                                   
%     x   :: A real input (exponent) value;                                                    
%     n   :: The number of terms of the MacLauring series to be included.                      
%                                                                                              
%  USES                                                                                        
%    double:: A built-in intrinsic function that converts an integer argument to a double precision value.  
%                                                                                              
%  REVISION DATE :: 03/18/2024                                                                 
%                                                                                              
%  ==================================================================================

term = 1.0;
sums = 1.0;

for k = 1:n-1
    term = term*x/double(k);
    sums = sums + term;
end

f = sums;
end