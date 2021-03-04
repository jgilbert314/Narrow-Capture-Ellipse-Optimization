function [ UB, LB ] = buildBounds(A, N)
% Author: Jason Gilbert
% Date: July 07, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function generates the upper and lower bounds of the optimization
%   space, in elliptical coordinates
% 
% Input: 
%   A - scaling factor determining the size of the elliptical domain
%   N - number of traps in configuration
% 
% Output:
%   LB - vector of lower bounds, in the form [distance, angle]
%        where distance = 0, and angle = 0
%   UB - vector of upper bounds, in the form [distance, angle]
%        where distance = A, and angle = 2*pi

vec = ones(1, N);
LB = [0*vec, 0*vec];
UB = [0.999*A*vec, 2*pi*vec]; % Max < A to keep traps from being placed on boundary

end