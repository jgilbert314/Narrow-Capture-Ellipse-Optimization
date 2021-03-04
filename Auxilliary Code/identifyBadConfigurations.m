function [ badN ] = identifyBadConfigurations(compArr, refArr, diff_thresh)
% Author: Jason Gilbert
% Date: February 22, 2020
% Last Updated: N/A
% Version: V00
% 
% Summary:
%   This function compares two sets of trap configurations and determines
%   which configurations need further improvement.
% 
% Input:
%   compArr, refArr - data to be compared to reference data: 
%       cell arrays of the format [N, x, p, t], where
%       N is the number of traps in each configuration
%       x is the coordinates of the traps
%       p is the interaction energy
%       t is the computation time
%       [1, 2:2*N+1, 2*N+2, 2*N+3] -> [N, x, p, t]
% 
% Output:
%   badN - an array of N corresponding to configurations which need further
%          optimization

% Extract N and p from logged data
refN = getN(refArr);
refP = getIntEng(refArr);
compN = getN(compArr);
compP = getIntEng(compArr);


% Extract N and p from
nComp = length(compArr);
badInd = false(1, nComp);
for itr = 1:nComp
    thisN = compN(itr);
    thisP = compP(itr);
    compN(itr) = thisN;
    
    % If the interaction energy > 0, mark config as bad. Otherwise, check
    % if any energys were less than the current.
    % Either means the configuration is bad. The latter is more time 
    % consuming to check.
    if (thisP > 0) % Optimal p is less than 0
        badInd(itr) = 1;
    else
        checkN = (refN <= thisN); % Indices of N <= thisN
        checkP = refP(checkN);    % p values corresponding to N <= thisN
        checkDiff = ( thisP - checkP );
        if ( any(checkDiff > diff_thresh) ) % Optimal p is minimal
            badInd(itr) = 1;
        elseif ( length(checkP) > 1 ) % p decreases as N increases
            if (any(checkP(1:end-1) < thisP))
                badInd(itr) = 1;
            end
        end
    end
end

badN = compN(badInd);

end
