function LinearizeRegress(ndata, x, y, model)
%  ==================================================================================
%  CODE7.2-LINEARIZE_REGRESS.m. A Matlab script module implementing Pseudocode 7.2.               
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
%  DESCRIPTION: A module to obtain least-squares-fit to the "power model" only. Users                 
%    can likewise incorporate the other linearizable models into the module.                   
%                                                                                              
%  ON ENTRY                                                                                    
%     n   :: The number of data in the set;                                                    
%    x,y  :: Arrays of length n containing the data;                                           
%   model :: Model flag, Model = 1 corresponds to the power model, Y= a0*x^b0.                 
%                                                                                              
%  ON EXIT                                                                                     
%    a0,b0:: Model parameters;                                                                 
%     E   :: Sum of the Squares of Residuals (SSR);                                            
%     S   :: Sum of the Squares of Mean Deviation (SSMD);                                      
%    r2   :: r-squared, coefficient of determination.                                          
%                                                                                              
%  USES                                                                                        
%     EXP :: Built-in Intrinsic function returning exponential of a real value, e^x.           
%     log :: Built-in Intrinsic function returning the natural log of a real value.            
%                                                                                              
%  REVISION DATE :: 03/03/2025                                                                 
%  ==================================================================================
    xx = zeros(1, ndata);
    yy = zeros(1, ndata);
    b = zeros(1, 2);
    c = zeros(2, 2);
    DD = 0; D1 = 0; D2 = 0; yavg = 0; yi = 0; xk = 0;

    for k = 1:ndata
        if model == 1
            xx(k) = log(x(k));
            yy(k) = log(y(k));
        else
            error('Undefined model...');
        end
    end

    for i = 1:2
        for j = 1:2
            p = i + j - 2;
            for k = 1:ndata
                xk = (p ~= 0) * xx(k)^p + (p == 0);
                c(i, j) = c(i, j) + xk;
            end
        end
        for k = 1:ndata
            p = i - 1;
            xk = (p ~= 0) * xx(k)^p + (p == 0);
            b(i) = b(i) + xk * yy(k);
        end
    end

    DD = c(1, 1) * c(2, 2) - c(2, 1) * c(1, 2);
    D1 = b(1) * c(2, 2) - b(2) * c(1, 2);
    D2 = b(2) * c(1, 1) - b(1) * c(2, 1);

    a0 = D1 / DD;
    b0 = D2 / DD;

    yavg = mean(y);

    S = 0;
    E = 0;

    if model <= 2
        a0 = exp(a0);
    else
        error('Undefined model...');
    end

    for k = 1:ndata
        S = S + (y(k) - yavg)^2;
        yi = a0 * x(k)^b0;
        E = E + (yi - y(k))^2;
    end

    r2 = 1 - E / S;

    fprintf('******* Best-Fit Coefficients and Parameters ********\n');
    fprintf('  a = %.6f   b = %.6f\n', a0, b0);
    fprintf('  E = %.6f   S = %.6f   r-squared = %.6f\n', E, S, r2);
end

% ==============================================================================
%  The main function to test LinearizeRegress.m
% ==============================================================================
ndata = 6;
x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
y = [5.3, 6.5, 6.8, 7.1, 7.8, 7.5];
model = 1;

LinearizeRegress(ndata, x, y, model);