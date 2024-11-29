%  ==================================================================================
%  CODE5.1-EXAMPLE5-4.M. A Matlab script module implementing Pseudocode 5.1.                      
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
%  DESCRIPTION: A program that calculates the first and second derivatives of                
%    a set of position (discrete) data.                                                                   
%                                                                                              
%  INPUT VARIABLES                                                                             
%     dt  :: Time increment (in seconds);                                                      
%     s   :: An array of length n, providing the distance traveled in meters.                  
%                                                                                              
%  OUTPUT VARIABLES                                                                            
%     v   :: Array of length n containing velocities (m/s) at discrete points;                 
%     a   :: Array of length n containing acelerations (m2/s) at discrete points.              
%                                                                                              
%  REVISION DATE :: 06/13/2024                                                                 
%  ================================================================================== 
clear; clc;

n = 6;
dt = 5.0;
s(1) = 0.0; s(2) = 5.45; s(3) = 21.3; s(4) = 82.84; 
s(5) = 212.86; s(6) = 473.6;

v(1) = (-s(3) + 4.0 * s(2) - 3.0 * s(1)) / (2.0 * dt);
a(1) = (-s(4) + 4.0 * s(3) - 5.0 * s(2) + 2.0 * s(1)) / (dt * dt);

% Apply central difference formulas
for i =2:n-1
    v(i) = (s(i + 1) - s(i - 1)) / (2.0 * dt);
    a(i) = (s(i + 1) - 2.0 * s(i) + s(i - 1)) / (dt * dt);
end

% Apply 2nd order backward difference formulas
v(n) = (s(n - 2) - 4.0 * s(n - 1) + 3.0 * s(n)) / (2.0 * dt);
a(n) = (-s(n - 3) + 4.0 * s(n - 2) - 5.0 * s(n - 1) + 2.0 * s(n)) / (dt * dt);

 % Print out the results
 t = 0.0;
 fprintf('Time (s)       s (m)       v (m/s)     a (m/s^2)\n');
 for i = 1:n
     fprintf("%9.4f   %9.4f   %9.4f   %9.4f\n", t, s(i), v(i), a(i));
     t = t + dt;
 end

 



