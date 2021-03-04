function [ bool_out ] = checkFileExists(parent_dir, file_name)
% Author: Jason Gilbert
% Date: February 22, 2020
% Last Updated: N/A
% Version: V00
% 
% Summary:
%   Checks if a file exists within a given directory.
% 
% Input: 
%   parent_dir - path to parent directory
%   file_name - name of file to look for in parent directory
% 
% Output:
%   bool_out - boolean variable indicating existance of file

ThisDir = dir(parent_dir); % Structure representation of directory
% If an item in ThisDir has the user-specified name, and is not a directory, 
% the directory exists
bool_out = any( strcmp({ThisDir.name}, file_name) & ~([ThisDir.isdir]) );

end