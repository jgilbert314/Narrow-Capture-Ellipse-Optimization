function [  ] = writeLogFile(dataArr, log_filepath)
% Author: Jason Gilbert
% Date: March 13, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function writes the contents of a cell array to file.
%   Each element of the cell array is written to a separate line, and the
%   file ends with a newline character.
%   If a file at the path specified already exists, it will be overwritten,
%   otherwise it will be created.
% 
% Input:
%   dataArr - a cell array which has vectors as elements, each of the format
%       [1, 2:2*N+1, 2*N+2, 2*N+3] -> [N, x, p, t]
%   log_filepath - path to where file should be written (including filename)
%       ex: <path>\<to>\<filename>.<ext>

nData = length(dataArr);

fileID = fopen(log_filepath, 'w');
try
    for itr = 1:nData
        
        thisN = dataArr{itr}(1); % Number of traps
        nvars = thisN*2;         % Number of coordinates
        thisCoord = dataArr{itr}(2:nvars+1); % Coordinates
        fval = dataArr{itr}(nvars+2);        % Interaction energy
        thisTime = dataArr{itr}(end);        % Computation time
        
        % Write to file
        fprintf(fileID, '%u,', thisN);
        fprintf(fileID, '%0.15f,', thisCoord);
        fprintf(fileID, '%0.15f,', fval);
        fprintf(fileID, '%0.15f\n', thisTime);
        
    end
catch ME
    fclose(fileID);
    rethrow(ME);
end
fclose(fileID);


end