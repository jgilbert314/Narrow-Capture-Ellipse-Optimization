function [A, nu] = calcCartesian2Elliptical(x, y, mu)
% Author: Jason Gilbert
% Date: July 21, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   Converts a set of cartesian coordinates to elliptical coordinates,
%   given a parameter which determines the eccentricity of the elliptical
%   domain, via the transformation
%       x = A*cosh(mu)*cos(nu)
%       y = A*sinh(mu)*sin(nu)
% 
% Input:
%   x - vector of x coordinates
%   y - vector of y coordinates
%   mu - constant determining the eccentricity of domain
% 
% Output:
%   A - vector of radial coordinates
%   nu - vector of angular coordinates

nu = atan2(y, x*tanh(mu)); % atan2 used to preserve sign
A = sqrt( (x.^2 + y.^2) ./ ((cosh(mu)*cos(nu)).^2 + (sinh(mu)*sin(nu)).^2) ); % This form used to avoid singularity at x=0 or y=0

end