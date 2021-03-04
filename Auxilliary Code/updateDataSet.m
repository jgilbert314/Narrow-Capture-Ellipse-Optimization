function [ writeData ] = updateDataSet(compData, refData)
% Author: Jason Gilbert
% Date: March 13, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function takes in two data sets and combines them into one.
%   If entries exist in both sets for the same N, the entry with the 
%   smaller interaction energy will be used. If an entry only exists for 
%   one of the two sets, it will be included in the output set.
% 
% Input:
%   compArr, refArr - two data sets to be combined into one
%       cell arrays of the format [N, x, p, t], where
%       N is the number of traps in each configuration
%       x is the coordinates of the traps
%       p is the interaction energy
%       t is the computation time
%       [1, 2:2*N+1, 2*N+2, 2*N+3] -> [N, x, p, t]
% 
% Output:
%   writeData - combined data set, in the same form as the input sets

% Extract relevant values from data sets
compN = getN(compData);
compP = getIntEng(compData);
refN = getN(refData);
refP = getIntEng(refData);

writeN = unique( [compN, refN] ); % Create list of all N
nWrite = length(writeN);
writeData = cell(1, nWrite);

for itr = 1:nWrite
   
    thisN = writeN(itr);
    compInd = find(compN == thisN);
    if ( isempty(compInd) ) % Entry unique to reference
        refInd = find(refN == thisN);
        writeData{itr} = refData{refInd};
    else
       refInd = find(refN == thisN);
       if ( isempty(refInd) ) % Entry unique to comparand
          writeData{itr} = compData{compInd}; 
       else % Entry found in both sets
          if ( compP(compInd) < refP(refInd) ) % Optimum in comparand set
              writeData{itr} = compData{compInd};
          else % Optimum in reference set
              writeData{itr} = refData{refInd};
          end
       end
    end
    
end

end