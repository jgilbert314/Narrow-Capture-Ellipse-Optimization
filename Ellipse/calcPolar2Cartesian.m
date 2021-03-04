function [x, y] = calcPolar2Cartesian(r, theta)
% Author: Jason Gilbert
% Date: December 07, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   Convertes vectors of radial and angular coordinates to vectors of 
%   Cartesian coordinates.
% 
% Input:
%   r - vector of radial coordinates
%   theta - vector of angular coordinates
% 
% Output:
%   x - horizontal coordinates
%   y - vertical coordinates

x = r.*cos(theta);
y = r.*sin(theta);

end