function [ AMFPT ] = calcAMFPT(omega, D, N, G_mat, nu, A)
% Author: Jason Gilbert
% Date: April 29, 2021
% Last Updated: N/A
% Version: V00
%
% Summary:
%   This function calculates the AMFPT for traps of common size using
%   equation (2.16), or for different sizes traps using equation (2.13).
%
% Input:
%   omega - area of domain: pi*a*b
%   D - diffusion coefficient
%   N - number of traps
%   G_mat - NxN matrix of pair-wise interaction terms. (eq 2.15)
%   nu - coefficients determined by size of traps -1/log(epsilon).
%        May be either a constant for traps of same size, or a vector for
%        traps of different size.
%   A - Nx1 vector of coefficients determined by trap locations
%
% Output:
%   AMFPT - the average mean-first passage timed


    if (size(nu) == 1) % Common trap size
        AMFPT = omega/(2*pi*D*nu*N) + 2*pi/N * ones(1, N)*G_mat*A;
    else % Different trap sizes
        AMFPT_vec = 2*pi*(G_mat*A) + A./nu;
        AMFPT = AMFPT_vec(1);
    end

end