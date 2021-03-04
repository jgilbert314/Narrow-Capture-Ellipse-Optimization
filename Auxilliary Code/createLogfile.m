function [ fileID, log_filename, log_filepath, fileInd ] = createLogfile(log_dir, log_filename_base, file_ext)
% Author: Jason Gilbert
% Date: July 07, 2020
% Version: V00
% Last Update: N/A
% 
% Summary: 
%   This function searches a given directory of logfiles in order to
%   specify an unused logfile name. A file with this name is then created 
%   and the information needed to identify the file in subsequent programs
%   is returned.
% 
% Input: 
%   log_dir - path to log files
%   log_filename_base - prefix string common to all log files
%   file_ext - string specifying log file extension
% 
% Output:
%   fileID - integer used to reference file in read/write operations
%   log_filename - name of log file
%   log_filepath - relative path to log file
%   fileInd - index specifying unique filename

fileInd = 0;
log_filename = [log_filename_base, num2str(fileInd), file_ext];
while (checkFileExists(log_dir, log_filename))
    fileInd = fileInd+1;
    log_filename = [log_filename_base, num2str(fileInd), file_ext];
end
log_filepath = [log_dir, '\', log_filename];

% Optimize Configurations
fileID = fopen(log_filepath, 'w'); % Create log file

end