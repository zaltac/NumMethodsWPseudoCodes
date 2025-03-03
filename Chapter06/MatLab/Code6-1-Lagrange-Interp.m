% ==============================================================================
%  The main program to test program LAGRANGE_EVAL.m
% ==============================================================================

n = 3;
x = [0.08, 0.25, 0.5, 0.90];
f = [0.25, 0.625, 0.81, 0.43];

xval = 0.1;
fval = LAGRANGE_EVAL(n, xval, x, f);
fprintf('fval = %f\n', fval);

function L = LAGRANGE_P(n, xval, x)
%  ==================================================================================
%  CODE6.1-LAGRANGE-P.m. A Matlab script module implementing Pseudocode 6.1.                    
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
%  DESCRIPTION: A module to compute the Lagrange polynomials at x=xval.                        
%                                                                                              
%  ON ENTRY                                                                                    
%    n    :: The number of data in the set minus 1;                                            
%    xval :: x value at which the dataset is to be interpolated;                               
%    x    :: Array of length n+1 containing abscissas, k=0,1,2,...,n.                          
%                                                                                              
%  ON EXIT                                                                                     
%    L    :: Array of length (n+1) containing Lagrange polynomials                             
%            that is, L(k) = L_k(xval) for k=0,1,2,...,n.                                      
%                                                                                              
%  REVISION DATE :: 02/28/2025                                                                 
%  ==================================================================================
L = ones(1, n+1);
for k = 0:n
    for i = 0:n
        if i ~= k
            L(k+1) = L(k+1) * (xval - x(i+1)) / (x(k+1) - x(i+1));
        end
    end
end
end

function fval = LAGRANGE_EVAL(n, xval, x, f)
%  ==================================================================================
%  CODE6.1-LAGRANGE_EVAL.m. A Matlab script module implementing Pseudocode 6.1.                      
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
%  DESCRIPTION: A function to evaluate Lagrange interpolating polynomial for an arbitrary xval 
%    within x0 <= xval <= xn; f=f(xval)=fval.                                                   
%                                                                                              
%  ARGUMENTS                                                                                   
%    n    :: The number of data in the set minus 1;                                            
%    xval :: x value at which the dataset is to be interpolated;                               
%    x    :: Array of length (n+1) containing abscissas, k=0,1,2,...,n;                        
%    f    :: Array of length (n+1) containing ordinates, k=0,1,2,...,n.                        
%                                                                                              
%  USES                                                                                        
%    LAGRANGE_P :: A subroutine generating the Lagrange polynomials.                                           
%                                                                                              
%  REVISION DATE :: 02/28/2025                                                                 
%  ==================================================================================

L = LAGRANGE_P(n, xval, x);

fval = sum(L .* f);
end
