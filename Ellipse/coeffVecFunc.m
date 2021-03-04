function [ A ] = coeffVecFunc(G_mat, N, omega, epsilon, D)
% Author: Jason Gilbert
% Date: December 25, 2019
% Last Updated: July 07, 2020
%   Updated to specifically handle the case epsilon = 0
% Version: V01
% 
% Summary:
%   This function calculates the vector of coefficients necessary to 
%   calculate the average mean first-passage time, as defined by 
%   equation (4.1).
% 
% Input:
%   G_mat - NxN matrix of pair-wise interaction terms. (eq 2.15)
%   N - number of traps
%   omega - area of domain: pi*a*b
%   epsilon - coefficient determining size of trap
%   D - diffusion coefficient
% 
% Output:
%   A - Nx1 vector of coefficients


I_mat = eye(N);             % Identity matrix

% Define right-hand side of equation
coeff = omega/2/pi/D/N;
RHS = ones(N, 1)*coeff;

if (epsilon > 0)
    nu = -1/log(epsilon);       % Trap-size dependent term
    E_mat = ones(N)/N;          % Square matrix of 1/N

    % Define operator acting on coefficient vector
    % left-hand side of equation
    LHS = ( I_mat + 2*pi*nu*(I_mat - E_mat)*G_mat );
else
    LHS = I_mat;
end

% Calculate vector of coefficients
A = LHS\RHS;

end