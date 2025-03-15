% ==============================================================================
%  The main program to test SUBROUTINE InverseL
% ==============================================================================
clear; clc;

n = 5;

L = [1.0 0.0 0.0 0.0 0.0; 
     4.0 3.0 0.0 0.0 0.0;
     5.0 2.0 4.0 0.0 0.0;
     3.0 1.0 8.0 2.0 0.0;
    -1.0 1.0 1.0 1.0 1.0];

fprintf(' Input Lower Matrix, L \n');
for i = 1:n
    fprintf('%.7f %.7f %.7f %.7f %.7f\n', L(i,1), L(i,2), L(i,3), L(i,4), L(i,5));
end

eL = InverseL(L, n);

fprintf(' ------ Output : Inverse L ------\n');
for i = 1:n
    fprintf('%.7f %.7f %.7f %.7f %.7f\n', eL(i,1), eL(i,2), eL(i,3), eL(i,4), eL(i,5));
end


function eL = InverseL(L, n)
%  ==================================================================================
%  CODE11.4-InverseL.m. A Matlab module implementing Pseudocode 11.4.                      
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
%   DESCRIPTION: A matlab module to invert a lower-triangular matrix.                                 
%                                                                                              
%   ON ENTRY                                                                                   
%     n    :: Dimension attribute of the matrix L;                                             
%     L    :: A lower-triangular matrix, nxn.                                                  
%                                                                                              
%   ON EXIT                                                                                    
%     eL   :: Inverse of L, also a lower-triangular matrix, nxn.                               
%                                                                                              
%   REVISION DATE :: 03/15/2025                                                                
%  ==================================================================================
    eL = zeros(n,n);
    
    eL(1,1) = 1/L(1,1);
    for i = 2:n
        eL(i,i) = 1/L(i,i);
        for j = i-1:-1:1
            sum = 0;
            for k = j+1:i
                sum = sum + eL(i,k)*L(k,j);
            end
            eL(i,j) = -sum/L(j,j);
        end
    end
end