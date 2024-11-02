% ==============================================================================
%  The main program to test mat_add.m
% ==============================================================================

% Initialize matrices A and B
A = [ 1  2  3 ;  4  5  6 ;  7  8  9];
B = [ 3 -2  1 ; -2  2  4 ;  3 -5  1];

% Print input matrix
fprintf('Input Matrix A');
disp(A)

fprintf('Input Matrix B\n');
disp(B)

% Perform matrix addition
fprintf('C = A + B Matrix\n');
C = mat_add(A, B);
disp(C)



function C = mat_add(A, B)  
%  ==================================================================================
%  CODE2.1-MAT_ADD.m. A Matlab script module implementing Pseudocode 2.1.                         
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
%   DESCRIPTION: A function to perform C=A+B matrix addition.                                
%                                                                                              
%   ON ENTRY                                                                                   
%     m,n :: Dimension attributes of the matrices;                                             
%      A  :: An input matrix (mxn);                                                            
%      B  :: An input matrix (mxn).                                                            
%                                                                                              
%   ON RETURN                                                                                    
%      C :: The output matrix (mxn).                                                           
%                                                                                              
%   REVISION DATE :: 03/18/2024                                                                
%  ==================================================================================     

mA = size(A, 1); nA = size(A, 2);
mB = size(B, 1); nB = size(B, 2);

if mA==mB && nA ==nB
    for i=1:mA
        for j=1:nA
            C(i,j) = A(i,j) + B(i,j);
        end
    end
else
    fprintf('A and B are not the same size!')
end
end


 