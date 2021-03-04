function [ p ] = meritFuncCircle(x_in)
% Author: Jason Gilbert
% Date: March 13, 2020
% Last Updated: N/A
% Version: V00
% 
% Summary:
%   This function calculates the interaction energy of a configuration of
%   traps.
% 
% Input:
%   x_in - a 1x2N vector of polar coordinates for N traps. The first N
%          elements correspond to radial coordinates, while the second N
%          correspond to angular coordinates
% 
% Output:
%   p - interaction energy


% Convert input to form used by greensMatFuncCircle
% Extract polar coordinates
N = length(x_in)/2;
% Convert to row vector
[xR, ~] = size(x_in);
if (xR ~= 1)
    x_in = x_in';
end
r = x_in(1:N);
theta = x_in(N+1:end);
% Convert to cartesian coordinates
x = r.*cos(theta);
y = r.*sin(theta);
% Form input vector
x_vec = [x ; y]';

    
% Calculate Green's function
G_mat = greensMatFuncCircle(x_vec, N);
% Calculate interaction energy
p = sum(sum(G_mat));


end