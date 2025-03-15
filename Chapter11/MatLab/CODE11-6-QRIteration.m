% ==============================================================================
%  The main program to test the module Basic_QR.m
% ==============================================================================
n = 4;
d = [ 3.0,  8.0, 6.0, 9.0];
e = [ 4.0,  2.0, 1.0, 0.0];
eps = 1.0e-3;
maxit = 299;
V = zeros(n, n);

[d, V] = Basic_QR(n, d, e, eps, maxit);

fprintf('\n*** EIGENVALUES ARE ***\n');
for i = 1:n
    fprintf('%d %f\n', i, d(i));
end
fprintf(' *** Eigenvectors *** \n');
for i = 1:n
    fprintf('%7.4f %7.4f %7.4f %7.4f\n', V(i,:));
end

function [d, V] = Basic_QR(n, d, e, eps, maxit)
%  ==================================================================================
%  CODE11.6-Basic_QR.m. A Matlab module implementing Pseudocode 11.6.                      
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
%   DESCRIPTION: A Matlab module implementing the QR Factorization algorithm to a symmetric       
%      tridiagonal matrix to find its eigenvalues and eigenvectors.                            
%                                                                                              
%   ON ENTRY                                                                                   
%      n   :: Dimension attribute of the tridiagonal matrix (nxn);                             
%      d   :: An array of length n containing the main diagonal, d(1) ... d(n);                
%      e   :: An array of length n containing the subdiagonal, e(1) ... e(n-1).                
%                                                                                              
%   ON EXIT                                                                                    
%      d   :: An array of length n containing the eigenvalues;                                 
%      V   :: A square matrix containing the eigenvector(nxn).                                 
%                                                                                              
%   USES                                                                                       
%     sqrt :: Built-in Intrinsic function returning the square root of a real value;           
%     zeros, eye are also built in functions.
%                                                                                              
%   REVISION DATE :: 03/15/2025                                                                
%  ==================================================================================
c = zeros(n-1, 1);
s = zeros(n-1, 1);

err = 1.0;
p = 0;
V = eye(n);

% ================ START THE QR ITERATION LOOP =====================
fprintf('*** Iteration History ***\n');
while (err > eps) && (p < maxit)
    t = e(1);
    for k = 1:n-1
        rho = sqrt(d(k)^2 + t^2);
        c(k) = d(k) / rho;
        s(k) = t / rho;
        d(k) = rho;
        t = e(k);
        e(k) = t*c(k) + d(k+1)*s(k);
        d(k+1) = -t*s(k) + d(k+1)*c(k);
        if k ~= n-1
            t = e(k+1);
            e(k+1) = t*c(k);
        end
        for i = 1:n
            q = V(i, k);
            r = V(i, k+1);
            V(i, k) = c(k)*q + s(k)*r;
            V(i, k+1) = -s(k)*q + c(k)*r;
        end
    end
    
    % *** Construct RQ Matrix   	   
    for k = 1:n-1
        d(k) = d(k)*c(k) + e(k)*s(k);
        t = d(k+1);
        e(k) = t*s(k);
        d(k+1) = t*c(k);
    end
    
    err = 0.0;
    for i = 1:n   % Calculate L2 norm of the vector of superdiagonal elements
        err = err + e(i)^2;
    end
    err = sqrt(err);
    p = p + 1;
    fprintf('iter=%3d %14.4e\n', p, err);
end

if p == maxit
    fprintf('Convergence within tolerance was not achived. Error is %14.4e\n', err);
end
% ********** END OF ITERATIONS ****************
end 