function [ G_mat ] = greensMatFunc(x_vec, a, b, f, omega, beta, xi_b)
% Author: Jason Gilbert
% Date: December 26, 2019
% Last Updated: N/A
% Version: V00
% 
% Input: 
%   x_vec - Nx2 array of trap coordinates (x, y)
%   a, b - length of major and minor axis of domain
%   f - domain dependent term: sqrt(a^2 + b^2)
%   omega - area of domain: pi*a*b
%   beta - domain dependent term: (a - b)/(a + b)
%   xi_b - a domain-dependent term: log( (a-b)/(a+b) )/2
% 
% Output:
%   G_mat - an NxN symmetric matrix of trap interaction terms defined by
%           equation (2.15)

[N, ~] = size(x_vec); % Number of traps
G_mat = zeros(N);  % Preallocate Green's Matrix


% Initialize terms independent of trap pairing
r_vec = sqrt( sum(x_vec.^2, 2) );      % Distance of traps from origin
xi_vec = zeros(1, N);
eta_vec = zeros(1, N);
for itr = 1:N
    x = x_vec(itr, 1);
    y = x_vec(itr, 2);
    
    xi_vec(itr) = xiFunc(x, y, f);
    eta_vec(itr) = etaFunc(x, y, f);

end


% Calculate terms dependent on trap pairing
for itrA = 1:N
    r = r_vec(itrA);
    xi = xi_vec(itrA);      % Equation (4.4a)
    eta = eta_vec(itrA);    % Equation (4.4b)
    for itrB = 1:itrA    
        r_0 = r_vec(itrB);
        xi_0 = xi_vec(itrB);
        eta_0 = eta_vec(itrB);
        
        z_arr = z_arrFunc(xi, xi_0, xi_b, eta, eta_0); % Equation (4.4b)
        
        % Greens function
        if (itrA ~= itrB)
            z_sum = z_sumFunc(beta, z_arr); % Term from Green's function
            xi_ang = max( [xi, xi_0] );
            G = greensFunc(r, r_0, a, b, omega, beta, xi_ang, z_sum);
            
            G_mat(itrA, itrB) = G;
            G_mat(itrB, itrA) = G;
        elseif (itrA == itrB)
            z_sum = zSing_sumFunc(beta, z_arr); % Term from Singular Green's Function
            G = greensSingFunc(r_0, a, b, omega, xi_0, eta_0, z_sum);
            G_mat(itrA, itrB) = G;
        end
                
    end
end


end