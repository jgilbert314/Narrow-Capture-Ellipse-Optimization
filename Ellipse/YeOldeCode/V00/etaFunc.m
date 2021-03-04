function [ eta ] = etaFunc(x, y, f)
% Author: Jason Gilbert
% Date: December 23, 2019
% Last Updated: N/A
% Version: V00
%
% Summary:
%   This function calculates one coordinate of the mapping of elliptical
%   domain to rectangular domain such that (x, y) -> (xi, eta),
%   defined by equation (4.4b)
% 
% Input:
%   x, y - point in elliptical coordinates
%   f - domain dependent term: sqrt(a^2 - b^2)
%
% Output:
%   eta - one coordinate of the mapping of elliptical domain to rectangular
%   domain such that (x, y) -> (xi, eta)


% Equation (4.4a)
mu = x^2 + y^2 - f^2;

% Equation (4.4b)
p = (-mu + sqrt(mu^2 + 4*f^2*y^2))/(2*f^2);
eta_st = asin(sqrt(p));

% Define eta accoring to quadrant
if ( (x >= 0) && (y >= 0) )     % Quadrant 1
    eta = eta_st;
elseif ( (x < 0) && (y >= 0) )  % Quadrant 2
    eta = pi - eta_st;
elseif ( (x <= 0) && (y < 0) )  % Quadrant 3
    eta = pi + eta_st;
elseif ( (x > 0) && (y < 0) )   % Quadrant 4
    eta = 2*pi - eta_st;
end

end