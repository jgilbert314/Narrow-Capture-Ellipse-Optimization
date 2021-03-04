function [ G ] = greensFunc(r, r_0, a, b, omega, beta, xi_ang, z_sum)
% Author: Jason Gilbert
% Date: December 23, 2019
% Last Updated: N/A
% Version: V00
%
% Summary:
%   This function calculates the Green's function, as defined by equation
%   (4.5a)
% 
% Input:
%   r, r_0 - distance of point (x, y) and (x_0, y_0) from origin
%   a, b - length of major and minor axis of domain
%   omega - area of domain: pi*a*b
%   beta - domain dependent term: (a - b)/(a + b)
%   xi_ang - max(xi, xi_0)
%   z_sum - approximation of the sum of the infinite series
% 
% Output:
%   G - value of Green's function, as defined by equation (4.5a)

% Green's Function
G = (r^2 + r_0^2)/4/omega - 3*(a^2 + b^2)/16/omega ...
    - log(beta)/4/pi - xi_ang/2/pi ...
    - z_sum/2/pi;

end