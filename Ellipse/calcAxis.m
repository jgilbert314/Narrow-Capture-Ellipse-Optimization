function [ a, b ] = calcAxis(ecc)
% Author: Jason Gilbert
% Date: June 02, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function calculates the lengths of the major and minor axis of an 
%   ellipse with a given eccentricity and an area of pi.
% 
% Input:
%   ecc - eccentricity of the ellipse, where 0 <= ecc <= 1
% 
% Output:
%   a - length of major axis, where a >= 1
%   b - length of minor axis, where b <= a

% Derived using the equations:
%   ecc = sqrt( 1 - (b/a)^2 )
%   A = pi*a*b
%       where A is the area of the ellipse, taken to be pi.

a = (1 - ecc^2)^-4; % Calc major axis
b = 1/a;            % Calc minor axis

end