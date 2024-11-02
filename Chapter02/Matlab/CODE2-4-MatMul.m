 
% PROGRAM Test_MAT_MUL
clear; clc;

% Define matrix dimensions
m = 2;
p = 4;
n = 3;

% Allocate memory for matrices

% Matrix A
A=[ 1.0 -3.0  5.0  1.0;  
   -2.0  4.0  1.0  2.0]

% Matrix B
B=[ 3.0   1.0  -2.0; 
    2.0   0.0   4.0;
   -1.0   1.0  -3.0;
    2.0   5.0   3.0 ]

fprintf('\n------ A(m,p)xB(p,n) matrix product is ------\n');
C = mat_mul(A, B);

fprintf('\nOUTPUT MATRIX C\n');
for i = 1:m
    fprintf('%.5f ', C(i, :));
    fprintf('\n');
end

 

function C = mat_mul(A, B)
%  ==================================================================================
%  CODE2.4-MAT_MUL.m. A Matlab script module implementing Pseudocode 2.4.                         
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
%   DESCRIPTION: A function to find A*B=C matrix multiplication.                             
%                                                                                              
%   ON ENTRY                                                                                   
%    m,p,n :: Dimension attributes of input/output matrices;                                   
%       A  :: An input matrix of size mxp;                                                    
%       B  :: An input matrix of size pxn.                                                    
%                                                                                              
%   ON RETURN                                                                                    
%       C  :: The output matrix of size mxn.                                                  
%                                                                                              
%  USES                                                                                        
%    size  :: Built-in function returning the size of an arrat;         
%    zeros :: Built-in function creating an mxn matrix of zeros.  
%
%   REVISION DATE :: 03/18/2024                                                                
%  ==================================================================================
    [m, pA] = size(A);
    [pB, n] = size(B);
    if pA==pB 
        C = zeros(m, n);  
	for i = 1:m
	    for j = 1:n
        	for k = 1:pA
	            C(i, j) = C(i, j) + A(i, k) * B(k, j);
	        end
	    end
	end
    else
       display("A and B are incompatible for matrix multiplication")
       display("Check your input arguments!..")
    endif
end