function [ xOpt, minDist ] = findDist2Bound(ecc, r0, x0)
% Author: Jason Gilbert
% Date: July 24, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   Given a point within an elliptical domain, with an area of pi, this 
%   function finds the minimum distance to the boundary from that point.
%   Minimum distance is found via fminsearch()
% 
% Input: 
%   ecc - eccentricity of elliptical domain
%   r0 - point in domain to compare to boundary, in the form [x, y]
%   x0 - initial point in x-domain to use for optimization search
% 
% Output:
%   xOpt - x-coordinate of closest point on boundary to given point
%   minDist - minimum distance to boundary from given point

[a, b] = calcAxis(ecc); % Calculate major/minor axis of ellipse with pi area

% Constrain problem to y >= 0, solution will be the same without having to
% check 2 cases (+/-)
if (r0(2) < 0)
    r0(2) = -r0(2);
end

% Function defining boundary, 
boundFunc = @(x) b*sqrt( 1 - (x/a).^2 );

% Function defining distance from arbitrary point r0 to boundary
distFunc = @(x) norm( [x,  boundFunc(x)] - r0 );

% Find x which minimizes distance to boundary
[xOpt, minDist] = fminsearch(distFunc, x0);


end