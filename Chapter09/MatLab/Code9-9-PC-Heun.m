% ==============================================================================
%  The main program to test PC_Heun.C
% ==============================================================================
x0 = 0.0; 
y0 = 2.0;
h = 0.1;
xlast = 1.0;

PC_Heun(h, x0, y0, xlast);


function PC_Heun(h, x0, y0, xlast)

%  ==================================================================================
%  CODE9.9-PC_HEUN.m. A Matlab script module implementing Pseudocode 9.9.                 
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
%  DESCRIPTION: A function module to estimate the solution of a first order IVP on [x0,xlast]           
%    using the Heun's Predictor-Correcter method. Numerical estimates are printed out, not stored.                       
%                                                                                              
%   ON ENTRY                                                                                   
%    h     :: Step size (it must be uniform);                                                  
%    x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps
%    xlast :: End point of the solution interval.                                              
%                                                                                              
%   Other Internal Variables                                                                   
%    x,y   :: Current estimates, x^(p+1) and y^(p+1).                                          
%                                                                                              
%  USES                                                                                        
%    abs  :: Built-in library function returning the absolute value of a real value.         
%    FCN  :: User-defined external function providing y'=f(x,y).                               
%                                                                                              
%  REVISION DATE :: 03/07/2025                                                                 
%  ==================================================================================
    x = x0;
    y = y0;
    
    fprintf('%12.7f %12.7f %14.3e\n', x, y, abs(y - exact(x)));
    
    while x < xlast
        % ==== PREDICTOR STEP 
        k1 = h * FCN(x, y);
        ys = y + k1;
        
        % ==== CORRECTOR STEP 
        x = x + h;
        k2 = h * FCN(x, ys);
        y = y0 + 0.5 * (k1 + k2);
        
        err = abs(y - exact(x));
        fprintf('%12.7f %12.7f %14.3e\n', x, y, err);
        
        x0 = x;
        y0 = y;
    end
end

function f = FCN(x, y)
% ==============================================================================
% DESCRIPTION: A function subprogram providing y'=f(x,y)
%
% ARGUMENTS:
%      x, y  :: Real input values.
%
% ==============================================================================
f = -y /(x + 1.0);
end

function y = exact(x)
% ==============================================================================
% DESCRIPTION: A function subprogram providing the true solution y=f(x) for 
%    testing the module. 
%
% ARGUMENTS:
%      x   :: A real input, independent variable.
%
% ==============================================================================
y = 2.0 / (x + 1.0 );
end
