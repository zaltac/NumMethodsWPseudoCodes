% ==============================================================================
%  The main program to test Richardson.M
% ==============================================================================
clear; clc;

x0 = 150.0; 
h = 5.0;
err = 1.0;
eps = 1e-6;

[D, nr, deriv] = Richardson(x0, h, eps);

for k = 1:nr
    fprintf('%2d %17.11g %17.11g %17.11g %17.11g %17.11g', k, D(k,1:k));
    fprintf('\n');
end

fprintf(' ---------------------------------\n');
fprintf('Derivative is= %19.11g\n', deriv);
fprintf(' ---------------------------------\n');

function [D, nr, deriv] = Richardson(x0, h, eps)
%  ==================================================================================
%  CODE5.2-RICHARDSON.M. A Matlab script module implementing Pseudocode 5.2.                      
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
%  DESCRIPTION: A function to compute the first derivative of an explicitly                      
%      defined function using Richardson's extrapolation.                                    
%                                                                                              
%  ON ENTRY                                                                                    
%      x0  :: Point at which derivative is to be computed;                                     
%      h   :: Initial interval size;                                                           
%      eps :: Tolerance desired.                                                               
%                                                                                              
%  ON RETURN                                                                                     
%      D   :: A matrix containing the Richardson's table (0..n, 0..n)                        
%      nr  :: Size of the table;                                                               
%    deriv :: Estimated derivative.                                                            
%                                                                                              
%  USES                                                                                        
%     ABS  :: Built-in Intrinsic function returning the absolute value of a real value.       
%                                                                                              
%  ALSO REQUIRED                                                                               
%     FUNC  :: User-defined external function providing the nonlinear equation.                
%                                                                                              
%  REVISION DATE :: 06/13/2024                                                                 
%  ==================================================================================

m = 1;
k = 1;
err = 1.0;
D = zeros(11,11);

while err > eps
    D(k,1) = (FUNC(x0+h) - FUNC(x0-h)) / (2.0*h);         % 1st derivative
    for m = 2:k
        D(k,m) = (4^(m-1) * D(k,m-1) - D(k-1,m-1)) / (4^(m-1) - 1);
    end
    if k >= 2  % Estimate diagonalwise differentiation error
        err = abs(D(k,k) - D(k-1,k-1));
    end
    h = h/2.0;
    k = k + 1;
end        
nr = k - 1;
deriv = D(nr,nr);
end

% ==============================================================================
%  USER-DEFINED FUNCTION "FUNC" OF ONE-VARIABLE
% ==============================================================================
function f = FUNC(x)
f = 25000/(-57.0 + x) - 5.2e6/x^2;
end

