% PROGRAM Test_Inv_Mat
clear; clc;
n = 5;
a = [1.0, 2.0, 1.0, 2.0, 3.0; 
     11.0, -1.0, 1.0, 4.0, 1.0;
      4.0, -1.0, 1.0, 1.0, -1.0;
     -3.0, 1.0, -8.0, -1.0, 5.0;
     -1.0, 1.0, 1.0, 1.0, 1.0];

for i = 1:n
    fprintf('%.5f   %.5f   %.5f   %.5f   %.5f\n', a(i,1), a(i,2), a(i,3), a(i,4), a(i,5));
end

b = a; % Save a for checking

fprintf('\n********* End of Input data *********\n\n');

ai = inv_mat(a);

fprintf('------ matrix A(-1) ---------------\n');
for i = 1:n
    fprintf('%.5f   %.5f   %.5f   %.5f   %.5f\n', ai(i,1), ai(i,2), ai(i,3), ai(i,4), ai(i,5));
end
fprintf('-----------------------------------\n');

c = b * ai;

fprintf('------ matrix A*A(-1) -------------\n');
for i = 1:n
    fprintf('%.5f   %.5f   %.5f   %.5f   %.5f\n', c(i,1), c(i,2), c(i,3), c(i,4), c(i,5));
end
fprintf('-----------------------------------\n');


%  ==================================================================================
%  CODE2.6-INV_MAT.m. A Matlab script module implementing Pseudocode 2.6.                         
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
%  DESCRIPTION: A function to find the inverse of a square matrix (with no pivoting).        
%                                                                                              
%  ON ENTRY                                                                                    
%     n  :: Dimension attribute of input matrix A;                                             
%     a  :: An input matrix (nxn).                                                             
%                                                                                              
%  ON RETURN                                                                                     
%     ai  :: Inverse of A(nxn).        
%                                                                                              
%  USES                                                                                        
%    size  :: Built-in function returning the size of an arrat;         
%    eye   :: Built-in function creating an identity matrix of size n.                                                         
%                                                                                              
%  REVISION DATE :: 03/18/2024                                                                 
%  ==================================================================================
function ai = inv_mat(a)
    n = size(a, 1);
    ai = eye(n);
    for j = 1:n
        p = 1/a(j,j);
        for k = 1:n
            a(j,k) = p*a(j,k);
            ai(j,k) = p*ai(j,k);
        end
        for i = 1:n
            s = a(i,j);
            if i ~= j
                for k = 1:n
                    a(i,k) = a(i,k) - s*a(j,k);
                    ai(i,k) = ai(i,k) - s*ai(j,k);
                end
            end
        end
    end
end


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
%   DESCRIPTION: A subroutine to find A*B=C matrix multiplication.                             
%                                                                                              
%   ON ENTRY                                                                                   
%    m,p,n :: Dimension attributes of input/output matrices;                                   
%       A  :: An input matrix of size mxp;                                                    
%       B  :: An input matrix of size pxn.                                                    
%                                                                                              
%   ON RETURN                                                                                    
%       C  :: The output matrix of size mxn.                                                  
%                                                                                              
%   REVISION DATE :: 03/18/2024                                                                
%  ==================================================================================
function c = mat_mul(a, b)
    [m, p] = size(a);
    [p, n] = size(b);
    c = zeros(m, n);
    for i = 1:m
        for j = 1:n
            for k = 1:p
                c(i,j) = c(i,j) + a(i,k) * b(k,j);
            end
        end
    end
end

