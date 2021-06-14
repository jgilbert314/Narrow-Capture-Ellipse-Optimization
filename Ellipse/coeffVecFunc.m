function [ A ] = coeffVecFunc(G_mat, N, omega, nu, D)
% Author: Jason Gilbert
% Date: December 25, 2019
% Last Updated: April 25, 2020
%   Updated to handle traps of different sizes. Function takes in nu
%   instead of epsilon, handles vector/singleton arguments for nu.
%   July 07, 2020
%   Updated to specifically handle the case epsilon = 0
% Version: V02
%
% Summary:
%   This function calculates the vector of coefficients necessary to
%   calculate the average mean first-passage time, as defined by
%   equation (4.1). 
%   In the case that traps are of the same size, equation (2.16) is used.
%   If the traps differ in size, equation (2.13) is used.
%   
%
% Input:
%   G_mat - NxN matrix of pair-wise interaction terms. (eq 2.15)
%   N - number of traps
%   omega - area of domain: pi*a*b
%   nu - coefficients determined by size of traps -1/log(epsilon)
%   D - diffusion coefficient
%
% Output:
%   A - Nx1 vector of coefficients


% Sum of elements of A is equal to this constant
coeff = omega/2/pi/D;


if (size(nu) == 1) % Common trap size
    
    % Define right-hand side of equation
    RHS = ones(N, 1)*coeff/N;
    
    I_mat = eye(N);             % Identity matrix
    
    if (nu ~= 0)
        % Define operator acting on coefficient vector
        % left-hand side of equation
        E_mat = ones(N)/N;          % Square matrix of 1/N
        
        LHS = ( I_mat + 2*pi*nu*(I_mat - E_mat)*G_mat );
    else
        LHS = eye(N); % Identity matrix
    end
    
    % Calculate vector of coefficients
    A = LHS\RHS;
    
    
else % Different trap sizes
    
    % Define right-hand side of equation
    RHS = nu;
    
    % Define left-hand side of equation
    LHS = zeros(N);
    for itrA = 1:N
        for itrB = 1:itrA
            
            term = 2*pi*nu(itrA)*G_mat(itrA, itrB);
            if (itrA == itrB)
                term = term + 1;
            end
            LHS(itrA, itrB) = term;
            LHS(itrB, itrA) = term; % Matrix is symmetric
            
        end
    end
    
    % Calculate vector of coefficients
    A_prop = LHS\RHS;                % Calculate A up to scaling factor
    corr_factor = sum(A_prop)/coeff; % Determine scaling factor
    A = A_prop/corr_factor;          % Correct for scaling factor
    
end


end