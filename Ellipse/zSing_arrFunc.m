function [ z_arr ] = zSing_arrFunc(xi, xi_0, xi_b, eta, eta_0)
% Author: Jason Gilbert
% Date: December 23, 2019
% Last Updated: N/A
% Version: V00
%
% Summary:
%   This function defines an array of values used to calculate the
%   Singular Green's function.
%   This array is defined in equation (4.6b)
%
% Input:
%   xi - mapping of ellipse to rectangle such that (x, y) -> (xi, eta)
%   xi_0 - mapping of reference point, in the same manner as xi
%   xi_b - a domain-dependent term: log( (a-b)/(a+b) )/2
%   eta - mapping of ellipse to rectangle such that (x, y) -> (xi, eta)
%   eta_0 - mapping of reference point, in the same manner as eta
%
% Output:
%   z_arr - a 1x8 array as defined by equation (4.6b)


% Intermediate terms
xi_p = xi + xi_0;
xi_m = abs( xi - xi_0 );
eta_p = 1i*( eta + eta_0 );
eta_m = 1i*( eta - eta_0 );

% Array elements
z_arr = zeros(1, 8); % Preallocate array
z_arr(1) = ( -xi_m + eta_m );
z_arr(2) = (  xi_m - 4*xi_b + eta_m );
z_arr(3) = ( -xi_p - 2*xi_b + eta_m );
z_arr(4) = (  xi_p - 2*xi_b + eta_m );
z_arr(5) = (  xi_p - 4*xi_b + eta_p );
z_arr(6) = ( -xi_p + eta_p );
z_arr(7) = (  xi_m - 2*xi_b + eta_p );
z_arr(8) = ( -xi_m - 2*xi_b + eta_p );
z_arr = exp(z_arr); % Vectorize exponentiation for efficiency

end