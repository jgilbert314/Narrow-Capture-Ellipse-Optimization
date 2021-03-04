function [ N ] = getN(dataArr)
% Author: Jason Gilbert
% Date: ??
% Version: V00
% Last Updated: N/A
% 
% Summary: 
% 	This function extracts the number of traps associated with each set
%   of configuration data. Data assumed to be in the form of a cell array,
%   with each element containing a vector in the form [N, ...]
% 
% Input:
%   dataArr - cell array with elements in the form of vectors [N, ...]
% 
% Output:
%   N - a vector of trap numbers, corresponding to each element of dataArr

nVal = length(dataArr);
N = zeros(1, nVal);
for itr = 1:nVal
    N(itr) = dataArr{itr}(1);
end


end
