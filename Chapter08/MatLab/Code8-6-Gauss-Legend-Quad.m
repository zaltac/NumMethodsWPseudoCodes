% ==============================================================================
%  The main program to test function Gauss_Legendre_Quad.m
% ==============================================================================
clear; clc;

n = 5;
x = zeros(n,1);
w = zeros(n,1);

eps = 1e-6;

[x, w] = Gauss_Legendre_Quad(n, eps);

fprintf('  i         x_i             w_i\n');
fprintf('--- ------------------- -------------------\n');
for i = 1:n
    fprintf('%2d %18.12f %18.12f\n', i, x(i), w(i));
end


%  ==================================================================================
%  CODE8.6-GAUSS_LEGENDRE_QUAD.m. A Matlab script module implementing Pseudocode 8.6.             
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
%  DESCRIPTION: A subroutine to generate N-point Gauss-Legendre quadrature                     
%    abscissas and weights on [-1,1].                                                          
%                                                                                              
%  ON ENTRY                                                                                    
%     n   :: Number of quadrature points;                                                      
%     eps :: Tolerance, i.e., desired level of numerical accuracy.                             
%                                                                                              
%  ON EXIT                                                                                     
%     x   :: Array of length N containing the abscissas;                                       
%     w   :: Array of length N containing the weights.                                         
%                                                                                              
%  USES                                                                                        
%    abs  :: Built-in Intrinsic function returning the absolute value of a real value;         
%    cos  :: Built-in Intrinsic function returning trig cosine value.  
%                                                                                              
%  REVISION DATE :: 03/04/2025                                                                 
%  ==================================================================================
function [x, w] = Gauss_Legendre_Quad(n, eps)
m = floor((n+1)/2);
pi = 3.1415926535897932385;

for i = 1:m
    u = cos(pi*(4*i-1)/(4*n+2));
    del = 1;
    while abs(del) > eps
        P0 = 1;
        P1 = u;
        for k = 2:n
            P2 = (2*k-1)*u*P1 - (k-1)*P0;
            P2 = P2/k;
            P0 = P1;
            P1 = P2;
        end
        PP = n*(u*P1 - P0)/(u^2 - 1);
        del = P1/PP;
        u = u - del;
    end
    x(i) = -u;
    w(i) = 2/(1-u^2)/PP^2;
    x(n+1-i) = u;
    w(n+1-i) = w(i);
end
end