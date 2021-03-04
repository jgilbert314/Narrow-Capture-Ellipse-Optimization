function [ x0 ] = buildInitialCoords(A, mu, N)
% Author: Jason Gilbert
% Date: July 07, 2020
% Version: V01
% Last Updated: December 07, 2020
%   Changed to handle case of circular domain
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
% 
% Output:
%   x0 - 1x2N vector in the form:
%           [r, nu]
%        where r are analagous to radial coordinates and nu are analagous
%        to angular coordinates.

r_0 = 0.1; % Min radius
r_M = 0.8; % Max radius

nu = mod(linspace(0, N, N), 2*pi); % Angular coordinate

if ( (A == 0) && isinf(mu) ) % If the domain is circular
    r = linspace(r_0, r_M, N);
else % If the domain is elliptical
    r = linspace(r_0*A, r_M*A, N).*sqrt( (cosh(mu)*cos(nu)).^2 + (sinh(mu)*sin(nu)).^2 ); % Wrong, but not very important
end
x0 = [ r, nu ];


end