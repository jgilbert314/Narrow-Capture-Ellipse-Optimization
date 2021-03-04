function [ file_path ] = formFilepath( folder_list, path_delim )
% Author: Jason Gilbert
% Date: August 12, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function creates a string representing a file path from a list of 
%   strings representing the folders in the path
% 
% Input: 
%   folder_list - cell array of strings to be joined
%   path_delim - character seperating directories in path
% 
% Output:
%   file_path - string representing file path

file_path = [];
for itr = 1:length(folder_list)-1
    file_path = [ file_path, folder_list{itr}, path_delim ];
end
file_path = [ file_path, folder_list{end} ];

end