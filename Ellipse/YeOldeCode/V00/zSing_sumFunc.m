function [ z_sum ] = zSing_sumFunc(beta, z_arr)
% Author: Jason Gilbert
% Date: June 02, 2020
% Last Updated: N/A
% Version: V00
%
% Summary:
%   This function approximates the sum of the infinite series appearing in 
%   the singular Green's function defined in equation (4.6a). This function
%   is distinct from the non-singular case in that the product of terms in
%   the array must be calculated excluding the first term, to avoid
%   evaluating log(0).
% 
% Input:
%   beta - domain dependent term: (a - b)/(a + b)
%   z_arr - array defined by equation (4.6b). This array may be obtained by
%           evaluating equation (4.5b).
% 
% Output:
%   z_sum - approximation of the sum of the infinite series appearing in 
%           equation (4.6a)

z_subarr = z_arr(2:end); % Subset of array

tol = 1e-15;    % Threshold determining convergence of infinite series
z_term = Inf;   % Initialize exit condition
z_sum = 0;      % Initialize output 
itr = 0;
while ( abs(z_term) > tol)
    
    beta_term = beta^(2*itr);
    z_term = log( prod( abs(1 - beta_term*z_subarr) ) );
    if (itr >= 1) % Component summed from n=1 to inf
       z_term = z_term + log(1 - beta_term);
    end
    z_sum = z_sum + z_term;
    
    itr = itr+1;
end


end