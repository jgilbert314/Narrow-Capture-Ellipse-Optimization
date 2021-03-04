function [ x0 ] = buildInitialCoords(A, mu, N)
% Author: Jason Gilbert
% Date: July 07, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function generates a set of elliptical coordinates to be used as
%   an initial point in the optimization
% 
% Input: 
%   A_max - scaling factor which determines the size of the elliptical
%           domain
%   mu - parameter which determines eccentricity
%   N - number of traps in configuration

nu = mod(linspace(0, N, N), 2*pi);
r = linspace(0.1*A, 0.8*A, N).*sqrt( (cosh(mu)*cos(nu)).^2 + (sinh(mu)*sin(nu)).^2 );
x0 = [ r, nu ];

end