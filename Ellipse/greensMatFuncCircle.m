function [ G_mat ] = greensMatFuncCircle(x_vec, N)
% Author: Jason Gilbert
% Date: ??, Documented December 07, 2020
% Last Updated: N/A
% Version: V00
% 
% Summary:
%   This function calculates a matrix of values produced by evaluating
%   the Green's function of a circle for all possible pairs of traps.
% 
% Input:
%   x_vec - Nx2 array of trap coordinates, where the first column contains
%           x-coordinates, and the second contains y-coordinates.
%   N - number of traps in the configuration
% 
% Output:
%   G_mat - NxN matrix of containing values of Green's function evaluated
%           for all possible trap pairs.

% Terms independent of trap pairing
r_vec = sqrt( sum(x_vec.^2, 2) ); % Distance from origin

G_mat = zeros(N);
for itrA = 1:N
    x = x_vec(itrA, :);
    r_x = r_vec(itrA);
    
    for itrB = 1:itrA
        xi = x_vec(itrB, :);
        r_xi = r_vec(itrB);
        
        % Terms appearing in Green's function
        r1 = sqrt(sum( (x - xi).^2 ));
        r2 = sqrt(sum( ( x*r_xi - xi/r_xi ).^2 ));
        r3 = r_x^2 + r_xi^2;
        
        % Calculate Green's function
        if (itrA ~= itrB)
            if (r_xi > 0)
                G = ( -log(r1) - log(r2) + r3/2 - 3/4 )/2/pi;
            else
                G = ( -log(r_x) + r_x^2/2 - 3/4 )/2/pi;
            end
            G_mat(itrA, itrB) = G;
            G_mat(itrB, itrA) = G; % Green's matrix is symmetric about diagonal
        else % Calculate singular Green's function
            if (r_xi > 0)
                c = sqrt(sum( (x*r_x - x/r_x).^2, 2 )); % Cross term
                R = ( -log(c) + r_xi^2 - 3/4 )/2/pi;
            else
                R = -3/8/pi;
            end
            G_mat(itrA, itrB) = R;
        end
        
    end
end


end