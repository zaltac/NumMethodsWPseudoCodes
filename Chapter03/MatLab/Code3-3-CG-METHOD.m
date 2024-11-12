% Test_CGM
% ==============================================================================
%  The main program to test FUNCTION CGM
% ==============================================================================

n = 10;
a = [2, -1,  0,  0,  0,  0,  0,  0,  0,  0;
    -1,  2, -1,  0,  0,  0,  0,  0,  0,  0;
     0, -1,  2, -1,  0,  0,  0,  0,  0,  0;
     0,  0, -1,  2, -1,  0,  0,  0,  0,  0;
     0,  0,  0, -1,  2, -1,  0,  0,  0,  0;
     0,  0,  0,  0, -1,  2, -1,  0,  0,  0;
     0,  0,  0,  0,  0, -1,  2, -1,  0,  0;
     0,  0,  0,  0,  0,  0, -1,  2, -1,  0;
     0,  0,  0,  0,  0,  0,  0, -1,  2, -1;
     0,  0,  0,  0,  0,  0,  0,  0, -1,  2];

b = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2.2]';

% Write input data
disp('    Coefficient Matrix    RHS    Initial Guess');

x = zeros(n, 1);  % SET INITIAL GUESS
for i = 1:n
    fprintf('%s %4.1f %4.1f\n', sprintf('%5.1f', a(i,:)), b(i), x(i));
end

eps = 1.0e-07;
maxit = 999;

% GO TO CGM MODULE
[x, iter, error] = CGM(n, eps, a, b, x, maxit);

% PRINT OUT THE RESULTS
disp(' ');
disp(' *** Solution ***');
for i = 1:n
    fprintf('%3d  %11.6f\n', i, x(i));
end

fprintf('\n Total no of iterations = %4d\n', iter);
fprintf(' Maximum Error          = %10.3e\n\n', error);


%  ==================================================================================
%  CODE3.3-CGM.M. A Matlab script module implementing Pseudocode 3.3.                             
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
%  DESCRIPTION: A subroutine to solve Ax=b linear system with the Conjugate Gradient Method.   
%                                                                                              
%  ON ENTRY                                                                                    
%     n   :: Number of equations;                                                              
%     A   :: Coefficient matrix (nÃ—n);                                                        
%     b   :: Array of length n containing the right hand side;                                 
%     x   :: Array of length n containing the initial guess;                                   
%    eps  :: Convergence tolerance;                                                            
%   maxit :: Maximum number of iterations.                                                     
%                                                                                              
%  ON RETURN                                                                                     
%     x   :: Array of length n containing the approximate solution;                            
%   iter  :: Total number of iterations performed;                                             
%   error :: Euclidean (L2-) norm of displacement at exit.                                     
%                                                                                              
%  USES                                                                                                                       
%    SQRT :: Built-in Intrinsic function returning the square root of a real value.                                         
%                                                                                              
%  REVISION DATE :: 12/11/2024                                                                 
%  ==================================================================================
function [x, iter, error] = CGM(n, eps, A, b, x, maxit)
    r = b - A * x;  % [r]^(0)=[b]-[A][x]^(0)
    d = r;          % [d]^(0)=[r]^(0)
    rho0 = r' * r;  % rho^(0)=[r^(0)].[r^(0)]
    
    disp(' ');
    disp('*** Iteration history ***');
    
    for p = 0:maxit
        Enorm = sqrt(rho0);  % E-norm=Sqrt([r]^(p-1).[r]^(p-1))
        fprintf('iter=%4d    E-norm= %11.5e\n', p, Enorm);
        
        if Enorm < eps
            break;  % CHECK FOR CONVERGENCE, EXIT if converged...
        end
        
        c = A * d;  % [c]^(p-1)=[A][d]^(p-1)
        rho = d' * c;  % rho=[d]^(p-1).[c]^(p-1)
        alpha = rho0 / rho;  % alpha^(p)=[r].[r]/([d].[c])
        x = x + alpha * d;  % [x]^(p)=[x]^(p-1)+alfa^(p)*[d]^(p-1)
        r = r - alpha * c;  % [r]^(p)=[r]^(p-1)-alfa^(p)*[d]^(p-1)
        rho = r' * r;  % rho^(p+1)=[r^(p)].[r^(p)]
        beta = rho / rho0;  % beta^(p) = rho/rho0
        d = r + beta * d;  % [d]^(p)=[r]^(p)+beta*[d]^(p-1)
        rho0 = rho;  % rho^(p)<==rho^(p+1)
    end
    
    iter = p;  % Set total number of iterations
    error = Enorm;  % Set recent Enorm as error
    
    if iter == maxit
        fprintf('\nFailed to converge after %d iterations\n\n', maxit);
    end
end