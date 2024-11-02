% PROGRAM Test_XdotY
clear; clc;
n = 4;
X = [1.0, 2.0, 3.0, 4.0];
Y = [4.0, 3.0, 2.0, 1.0];

fprintf('\n\t Input Vector X\n');
for i = 1:n
    fprintf('x(%d) = %f\n', i, X(i));
end

fprintf('\n\t Input Vector Y\n');
for i = 1:n
    fprintf('y(%d) = %f\n', i, Y(i));
end

fprintf('\n------ X*Y dot product is ==> %f\n', XdotY(n, X, Y));

function result = XdotY(n, x, y)
%  ==================================================================================
%  CODE2.3-XdotY.m. A Matlab script module implementing Pseudocode 2.3.                           
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
%   DESCRIPTION: A function to compute the dot product of two vectors, x and y.                
%                                                                                              
%   ARGUMENTS                                                                                  
%      n   :: Dimension attribute of the input vectors;                                        
%     x, y :: The input vectors of length n.                                                   
%                                                                                              
%   REVISION DATE :: 03/18/2024                                                                
%  ==================================================================================
sums = 0.0;
  for i = 1:n
     sums = sums + x(i)*y(i);
  end
  result = sums;
end