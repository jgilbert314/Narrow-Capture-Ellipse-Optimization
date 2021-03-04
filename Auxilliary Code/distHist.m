function [ hist2D, bins ] = distHist(res, r, chordList)
% Author: Jason Gilbert
% Date: April 25, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function takes a vector of radial coordinates, and a cell array
%   of distances from each trap to its neighbours, and produces a 2D 
%   histogram which sorts the traps according to their radial coordinates
%   (rows), and their distance from the reference trap (columns)
% 
% Inputs:
%   res - number of bins to use
%   r - a vector of trap radii
%   chordList - a 1xN cell array, where each element contains a vector of 
%               distances between a trap and its neighbours
% 
% Outputs:
%   Outputs are produced by the Matlab function histcounts().
%   See that documentation for more details.
%   hist - an MxM histogram, where M = res-1. Each row of the array 
%          contains the number of traps within the corresponding radial 
%          interval, while each column contains the number of traps within 
%          the corresponding interval of distances between traps.
%   bins - a vector of bin edges, defining the bin intervals via:
%               b_n = bins(n) - bins(n-1), where bins(0) = 0 ;
%          each of which has a width of 1/(res - 1)


% Assign radial coordinates to histogram bins defined by 'edges'. Return
% the indices of the bin each coordinate was assigned to, so that
% histograms of pair-wise distances can be assigned to the corresponding to
% trap coordinate
edges = linspace(0, 1, res+1); % Number of bins is one less than edges
bins = edges(2:end);
[~, ~, histInd] = histcounts(r, edges); 


% Build 2D histogram as a function of trap radii and distance to nearby
% traps
hist2D = zeros(res);
for itr = 1:length(histInd)
    thisInd = histInd(itr);
    thisChord = chordList{itr};
    
    % Sort traps according to their distance from reference trap. Add
    % result to corresponding radial bin
    thisHist = histcounts(thisChord, edges);
    hist2D(thisInd, :) = hist2D(thisInd, :) + thisHist;
end


end