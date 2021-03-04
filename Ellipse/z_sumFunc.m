function [ z_sum ] = z_sumFunc(beta, z_arr)
% Author: Jason Gilbert
% Date: December 23, 2019
% Last Updated: N/A
% Version: V00
%
% Summary:
%   This function approximates the sum of the infinite series appearing in 
%   the Green's function defined in equation (4.5a)
% 
% Input:
%   beta - domain dependent term: (a - b)/(a + b)
%   z_arr - array defined by equation (4.5b)
% 
% Output:
%   z_sum - approximation of the sum of the infinite series appearing in 
%           equations (4.5a) and (4.6a)

tol = 1e-15;    % Threshold determining convergence of infinite series
z_term = Inf;   % Initialize exit condition
z_sum = 0;      % Initialize output 
itr = 0;
while ( abs(z_term) > tol)
    
    z_term = log( prod( abs(1 - beta^(2*itr)*z_arr) ) );
    z_sum = z_sum + z_term;
    
    itr = itr+1;
end

