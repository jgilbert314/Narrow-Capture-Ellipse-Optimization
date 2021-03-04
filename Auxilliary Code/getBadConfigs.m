function [ badN ] = getBadConfigs(compN, compP, refN, refP, diff_thresh)
% Author: Jason Gilbert
% Date: February 23, 2020
% Last Updated: N/A
% Version: V00
% 
% Summary:
%   This function determines which trap configurations must be improved
%   
% Input:
%   compN - vector of trap numbers for configurations to be improved
%   compP - vector of interaction energies to be improved
%   refN - vector of trap numbers for reference configureations
%   refP - vector of interaction energies to be used as reference
%   diff_thresh - maximum ammount compared values may be greater than 
%                 reference
% 
% Output:
%   badN - numbers of trap configurations to be improved


% Extract N and p from
nComp = length(compN);
badInd = false(1, nComp);
for itr = 1:nComp
    thisN = compN(itr);
    thisP = compP(itr);
    compN(itr) = thisN;
    
    % If the interaction energy > 0, mark config as bad. Otherwise, check
    % if any energys were less than the current.
    % Either means the configuration is bad. The latter is more time 
    % consuming to check.
    if (thisP > 0)
        badInd(itr) = 1;
    else
        checkN = (refN <= thisN);
        checkDiff = ( thisP - refP(checkN) );
        if ( ( any( checkDiff > diff_thresh ) || (thisP > 0) ) )
            badInd(itr) = 1;
        end
    end
end

badN = compN(badInd);

end