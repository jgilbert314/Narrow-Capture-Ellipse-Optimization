function [ p ] = getIntEng(dataArr)
% Author: Jason Gilbert
% Date: ??
% Version: V00
% Last Updated: N/A
% 
% Summary: 
% 	This function extracts the number of traps associated with each set
%   of configuration data. Data is assumed to be in the form of a cell 
%   array, with vectors in each element
% 
% Input:
%   dataArr - cell array with elements in the form of vectors 
%             [1, ..., 2*N+2, ...] -> [N, ..., p, ...]
% 
% Output:
%   p - a vector of energies, corresponding to each element of dataArr

nVal = length(dataArr);
p = zeros(1, nVal);
for itr = 1:nVal
    thisN = dataArr{itr}(1);
    p(itr) = dataArr{itr}(2*thisN + 2);
end


end
