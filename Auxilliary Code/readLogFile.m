function [ dataArr ] = readLogFile(log_filename)
% Author: Jason Gilbert
% Date: ?? 
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function imports a .csv file as a a cell array, where each element
%   of the array is a vector which contains the elements of the 
%   corresponding row of text.
% 
% Input:
%   log_filename - name of file to be imported, including file extension.
%                  File should be comma-delimited and contain only
%                  numerical values.
% 
% Output:
%   dataArr - a 1xN cell array, where N is the number of rows in the file.
%             Each element of the array is a row vector containing the
%             elements of the corresponding row from the file.

% TODO: add optional argument for specifying delimiter.


% Count rows for preallocation
fileID = fopen(log_filename, 'r');
try
    rowCount = -1; % Start at -1 so the min value is 0
    line = 0;
    while (line ~= -1)
        line = fgets(fileID);
        rowCount = rowCount+1;
    end
    
catch ME
    fclose(fileID);
    rethrow(ME);
end
fclose(fileID);

% Import file text
data = cell(1, rowCount);
fileID = fopen(log_filename, 'r');
try
    for itr = 1:rowCount
        data{itr} = fgets(fileID);
    end
catch ME
    fclose(fileID);
    rethrow(ME);
end
fclose(fileID);



%% Convert strings to arrays
delim = ',';

dataArr = cell(1, rowCount);
for itr = 1:rowCount
    thisLine = data{itr};
    sepInds = find(thisLine == delim); % Find delimiter positions
    nSep = length(sepInds);            % Number of values
    
    thisArr = zeros(1, nSep+1); % Pre-allocate data
    thisArr(1) = str2double( thisLine(1:sepInds(1)) ); % Assign first value
    for itrR = 1:nSep-1 % Loop through body
        thisInd = sepInds(itrR);
        nextInd = sepInds(itrR+1);
        thisArr(itrR+1) = str2double( thisLine(thisInd:nextInd) );
    end
    thisArr(end) = str2double( thisLine(sepInds(end):end) ); % Assign last value
    
    dataArr{itr} = thisArr; % Assign array to returned cell array
    
end


end