function [ xi ] = xiFunc(x, y, f)
% Author: Jason Gilbert
% Date: December 23, 2019
% Last Updated: N/A
% Version: V00
%
% Summary:
%   This function calculates one coordinate of the mapping of elliptical
%   domain to rectangular domain such that (x, y) -> (xi, eta),
%   defined by equation (4.4a)
%
% Input:
%   x, y - point in elliptical coordinates
%   f - domain dependent term: sqrt(a^2 - b^2)
%
% Output:
%   xi - one coordinate of the mapping of elliptical domain to rectangular
%   domain such that (x, y) -> (xi, eta)

mu = x^2 + y^2 - f^2;
s = (-mu - sqrt( mu^2 + 4*f^2*y^2 ))/(2*f^2);
xi = log(1 - 2*s + 2*sqrt(s^2 - s))/2;

end