% PROGRAM Test_AX
n = 4;
i = 1:n;
j = 1:n;
A = [1.0, 2.0, 3.0, -2.0; 
     4.0, 1.0, 2.0,  3.0;
     3.0, 2.0, 1.0,  2.0;
    -2.0, 3.0, 4.0,  1.0];
x = [2.0, 5.0, 3.0, -3.0];

fprintf('\n    INPUT MATRIX A\n');
fprintf('%10.5f %10.5f %10.5f %10.5f\n', A');

fprintf('\n    INPUT VECTOR x\n');
for i = 1:n
    fprintf('x(%d) = %10.5f\n', i, x(i));
end

b = zeros(n, 1);
b = Ax(n, A, x);

fprintf('\n------ A(n,n)X(n) product is \n');

fprintf('\n    OUTPUT VECTOR b\n');
fprintf('%10.5f\n', b);

function b = Ax(n, A, x)
%  ==================================================================================
%  CODE2.5-Ax.m. A Matlab script module implementing Pseudocode 2.5.                              
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
%  DESCRIPTION: A function to perform A * x = b matrix-vector multiplication.                
%                                                                                              
%  ON ENTRY                                                                                    
%     n   :: Dimension attributes of input/output matrices;                                    
%     A   :: An input matrix of size nxn;                                                     
%     x   :: An input vector of length n.                                                      
%                                                                                              
%  ON EXIT                                                                                     
%     b   :: The output vector of length n.                                                    
%                                                                                              
%  REVISION DATE :: 03/18/2024                                                                 
%  ==================================================================================
for i = 1:n
   b(i) = 0.0;
   for j = 1:n
      b(i) = b(i) + A(i,j) * x(j);
   end
end

end
