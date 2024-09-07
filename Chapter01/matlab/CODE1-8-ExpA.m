% PROGRAM Test_EXPA
% ==============================================================================
%  The main program to test function EXPA
% ==============================================================================
x = input("Enter x: ");
eps = 1e-6;
fprintf('e^%.3f = %.7f\n', x, EXPA(x, eps));


function f = EXPA(x, eps)
%  ==================================================================================
%  CODE1.8-ExpA.m. A Matlab script module implementing Pseudocode 1.8.                            
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
%  DESCRIPTION: A function to compute e^x adaptively using the MacLaurin series                
%     within a user-defined tolerance.                                                         
%                                                                                              
%  ARGUMENTS                                                                                   
%     x   :: A real input (exponent) value;                                                    
%    eps  :: A user-defined convergence tolerance.                                             
%                                                                                              
%  USES                                                                                        
%   Float :: A built-in intrinsic function that converts an integer argument to a real value.  
%    Abs  :: Built-in Intrinsic function returning the absolute value of a real value.         
%                                                                                              
%  REVISION DATE :: 04/11/2024                                                                                                                                       
%  ================================================================================== 
k = 0;
term = 1.0;
sums = 1.0;
while abs(term) > eps * abs(sums)
    k = k + 1;
    term = term * x / k;
    sums = sums + term;
end
f = sums;
end