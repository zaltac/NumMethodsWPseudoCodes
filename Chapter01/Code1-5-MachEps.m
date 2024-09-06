% ==============================================================================
%  The main program to test FUNCTION MACHEPS
% ==============================================================================
%x = input("Enter a real value in 0<x<=1: ");
x=0.5;
% Print result
fprintf("Machine Epsilon is= %e\n", MACHEPS(x));
 

function epsm = MACHEPS(x)
%  ==================================================================================
%  CODE1.5-MachEps.m. A matlab script module for implementing Pseudocode 1.5.           
% 
%  NUMERICAL METHODS FOR SCIENTISTS AND ENGINEERS: WITH PSEUDOCODES
%  First Edition. (c) By Zekeriya ALTAÇ (2024).
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
%  DESCRIPTION: A function module to calculate the machine epsilon, which            
%     provides a positive machine value that is almost negligible compared to 1.     
%                                                                                    
%  INPUT ARGUMENT                                                                    
%     x   :: A real number, 0<x<=1.                                                  
%                                                                                    
%  REVISION DATE :: 03/18/2024                                                       
%                                                                                    
%  ==================================================================================
while (1.0 + x/2.0) > 1.0
    x = x/2.0;
end
epsm = x;
end
