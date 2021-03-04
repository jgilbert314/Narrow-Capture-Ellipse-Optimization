function [ R ] = greensSingFunc(r_0, a, b, omega, xi_0, eta_0, z_sum)
% Author: Jason Gilbert
% Date: December 23, 2019
% Last Updated: N/A
% Version: V00
%
% Summary:
%   This function calculates the singular part of the Green's function, as
%   defined by equation (4.6a)
% 
% Input:
%   r_0 - distance of point (x_0, y_0) from origin
%   a, b - length of major and minor axis of domain
%   omega - area of domain: pi*a*b
%   xi_0 - mapping of ellipse to rectangle such that (x_0, y_0) -> (xi_0, eta_0)
%   eta_0 - mapping of ellipse to rectangle such that (x_0, y_0) -> (xi_0, eta_0)
%   z_sum - approximation of the sum of the infinite series
% 
% Output:
%   R - value of Singular Green's function, as defined by equation (4.6a)

% TODO: Factor Green's function to minimize operations
%       Group terms by common factor (2*pi, omega, etc)

% Singular Green's Function
R = r_0^2/2/omega - 3*(a^2 + b^2)/16/omega + log(a + b)/2/pi ...
    - xi_0/2/pi + log( cosh(xi_0)^2 - cos(eta_0)^2 )/4/pi ...
    - z_sum/2/pi;

end