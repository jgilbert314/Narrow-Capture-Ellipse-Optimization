function [ eps_p ] = calcTwoTrapSizes(eps_0, N_0, N_p, k)
% Author: Jason Gilbert
% Date: April 29, 2021
% Last Updated: N/A
% Version: V00
%
% Summary:
%   For a set of traps of two different sizes, this function calculates
%   the size of traps needed to preserve the area of a reference absorbing
%   set, when the size of one trap is specified relative to the other.
%
% Input:
%   eps_0 - size of traps in reference set
%   N_0 - number of traps in reference set
%   N_p - number of traps differing in size
%   k - factor determining size of one trap set, relative to reference
%
% Output:
%   eps_p - 1xN vector defining trap sizes

% Calculated using the conditions that:
%   A_n = N_n*pi*eps_n^2
%   A_0 = A_1 + A_2
%   eps_1 = k*eps_0



if (N_0 == N_p) % Force k = 1
   eps_p = eps_0*ones(N_0, 1);
else
    eps_1 = k*eps_0;
    eps_2 = sqrt( (1 - k^2)*N_0/N_p + k^2 )*eps_0;
    
    eps_p = zeros(N_0, 1);
    N_1 = N_0 - N_p; % Number of traps of size eps_1
    eps_p(1:N_1) = eps_1;
    eps_p(N_1+1:end) = eps_2;
end

end