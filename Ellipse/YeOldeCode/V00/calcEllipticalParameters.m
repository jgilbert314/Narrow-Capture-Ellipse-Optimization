function [A, mu] = calcEllipticalParameters(a, b)
% Author: Jason Gilbert
% Date: July 07, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function calculates the paramters which define an ellipse in
%   elliptical coordinates given the major and minor axis of the ellipse,
%   where the elliptical coordinates are related to the cartesian
%   coordinates via:
%       x = A*cosh(mu)*cos(nu)
%       y = A*sinh(mu)*sin(nu)
%       where mu > 0, and 0 <= nu < 2*pi
% 
% Input:
%   a - major axis of ellipse
%   b - minor axis of ellipse
% 
% Output:
%   A - scaling constant determining size of ellipse
%   mu - constant determining the eccentricity of ellipse

% Found by evaluating at:
%   nu = 0, x = a, y = 0
%   nu = pi/2, x = 0, y = b
% and taking the ratio of the two equations
mu = atanh(b/a); 
% Found by evaluating at:
%   nu = pi/2, x = 0, y = b
A = b/sinh(mu);

end