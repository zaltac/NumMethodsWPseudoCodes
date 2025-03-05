% ==============================================================================
%  The main program to test Runge_Kutta.m
% ==============================================================================
clear; clc;
x0 = 0.0; 
xlast = 1;
y0 = 2.0;
h = 0.1;
n = input('Enter Order of RK method ');
Runge_Kutta(n, h, x0, y0, xlast);



function Runge_Kutta(n, h, x0, y0, xlast)
%  ==================================================================================
%  CODE9.4-Runge_Kutta.m. A Matlab script module implementing Pseudocode 9.4.                     
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
%    using the Runge-Kutta methods. Numerical estimates are printed out, not stored.           
%                                                                                              
%  ON ENTRY                                                                                    
%    n     :: Order of the Runge-Kutta scheme;                                                 
%    h     :: Step size (it must be uniform);                                                  
%    x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps
%   xlast  :: End point of the solution interval.                                              
%                                                                                              
%  Other Internal Variables                                                                    
%    x,y   :: Current estimates, x^(p+1) and y^(p+1).                                          
%                                                                                              
%   USES                                                                                       
%     abs  :: Built-in Intrinsic function returning the absolute value of a real value.        
%    DRV_RK:: A driver subprogram performing one-step RK scheme.                               
%                                                                                              
%  REVISION DATE :: 03/05/2025                                                                 
%  ==================================================================================
    fprintf('%12.7f %12.7f\n', x0, y0);
    
    x = x0;
    while x < xlast
        [x, y] = DRV_RK(n, h, x0, y0);
        yt = exact(x);
        aerr = abs(y - yt);
        fprintf('%12.7f %12.7f %12.7f %14.3e\n', x, yt, y, aerr);
        x0 = x; y0 = y;
    end
end

function [x, y] = DRV_RK(n, h, x0, y0)
%  ==================================================================================
%  CODE9.4-DRV_RK.m. A Matlab script module implementing Pseudocode 9.4.                          
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
%  DESCRIPTION: A driver function module employing one-step RK2, RK3, or RK4 scheme.               
%                                                                                              
%  ON ENTRY                                                                                    
%   n     :: Order of Runge-Kutta scheme;                                                      
%   h     :: Step size (it must be uniform);                                                   
%   x0, y0:: Initial values, also denotes prior estimates, x^(p) and y^(p), on following steps;
%                                                                                              
%  ON EXIT                                                                                     
%   x,y   :: Current estimates, x^(p+1) and y^(p+1).                                           
%                                                                                              
%  USES                                                                                        
%    ABS  :: Built-in Intrinsic function returning the absolute value of a real value.         
%    FCN  :: User-defined external function providing y'=f(x,y).                               
%                                                                                              
%  REVISION DATE :: 03/05/2025                                                                 
%  ==================================================================================
    hlf = 0.50;
    xh = x0 + 0.5 * h;
    x1 = x0 + h;
    
    switch n
        case 2 % Case of RK2
            xk1 = h * fcn(x0, y0);
            ym = y0 + xk1;
            xk2 = h * fcn(x1, ym);
            xk = hlf * (xk1 + xk2);
        case 3 % Case of RK3
            xk1 = h * fcn(x0, y0);
            ym = y0 + hlf * xk1;
            xk2 = h * fcn(xh, ym);
            ym = y0 - xk1 + 2 * xk2;
            xk3 = h * fcn(x1, ym);
            xk = (xk1 + 4 * xk2 + xk3) / 6;
        case 4 % Case of RK4
            xk1 = h * fcn(x0, y0);
            ym = y0 + hlf * xk1;
            xk2 = h * fcn(xh, ym);
            ym = y0 + hlf * xk2;
            xk3 = h * fcn(xh, ym);
            ym = y0 + xk3;
            xk4 = h * fcn(x1, ym);
            xk = (xk1 + 2 * xk2 + 2 * xk3 + xk4) / 6;
        otherwise
            fprintf(' PROGRAM DOES NOT HANDLE CASE OF N=%d\n', n);
            return;
    end
    y = y0 + xk;
    x = x1;
end

function f = fcn(x, y)
% ==============================================================================
% DESCRIPTION: A function module providing y'=f(x,y)
%
% ARGUMENTS:
%      x, y  :: Real input values.
%
% ==============================================================================
    f = -y / ( x + 1 );
end

function y = exact(x)
% ==============================================================================
% DESCRIPTION: A function module providing the true solution y=(x) for 
%    testing the module. 
%
% ARGUMENTS:
%      x   :: A real input, independent variable.
%
% ==============================================================================
    y = 2 / ( x + 1 );
end