function [x, y] = calcElliptical2Cartesian(A, mu, nu)
% Summary:
%   Converts a set of elliptical coordinates to cartesian coordinates via
%   the transformations:
%       x = A*cosh(mu)*cos(nu)
%       y = A*sinh(mu)*sin(nu)
%       where mu > 0, and 0 <= nu < 2*pi
% 
% Input:
%   A - vector of radial coordinates
%   mu - constant which determines eccentricity of ellipse
%   nu - vector of angular coordinate
% 
% Output:
%   x - vector of horizontal coordinates
%   y - vector of vertical coordinates

x = cosh(mu)*A.*cos(nu);
y = sinh(mu)*A.*sin(nu);

end