function [  ] = plotTrapDist(DataSet, ecc, markSize)
% Author: Jason Gilbert
% Date: July 21, 2020
% Version: V01
% Last Updated: December 20, 2020
%   - Added argument (markSize) for specifying marker sizes
% 
%   July 24, 2020
%   - Now properly calculates shortest distance between a trap and the
%   boundary
% 
% Summary: 
%   This function calculates the minimum, and mean, distances between
%   traps, as well as the minimum distance from a trap to the boundary, and
%   plot the results vs the number of traps in the domain.
% 
% Input:
%   DataSet - a structure of trap configuration data with the fields:
%               ConfigData  - a structure of objects which describe the
%                             trap configuration, with the following fields
%                   thisN - an integer defining the number of traps in the
%                           configuration
%                   coordVec - an array of trap coordinates, in the format:
%                              [x1, y1 ; ... ; xN, yN]
%                   chordList - a 1xM cell array containing vectors of
%                               distances between traps
%                   triD - an array of indices produced by delaunay() which
%                          defines the nearest neighbours of each trap
%               structIndex - an Mx2 array which defines a mapping of
%                             structure indices to trap numbers
%	ecc - eccentricity of domain. If unspecified, the
%         domain will be considered circular
%   markSize - size of markers used for plot
% 
% Output:
%   N/A

% Initialization
[a, b] = calcAxis(ecc); % Calculate the major/minor axis of the domain

numSets = length(DataSet.ConfigData); % Total number of data sets
distPerTrap = zeros(1, numSets);      % Diameter of circle occupying average area per trap
meanPair = zeros(1, numSets);         % Mean distance between traps
minPair = zeros(1, numSets);          % Minimum distance between traps
minBound = zeros(1, numSets);         % Minimum distance to boundary
N_vals = zeros(1, numSets);           % Number of traps per configuration


% Calculate distances
for itr = 3:numSets
    ThisConfig = DataSet.ConfigData(itr);
    N = ThisConfig.thisN;
    N_vals(itr) = N;
    
    % Calculate average area
    distPerTrap(itr) = 2/sqrt(N);
    
    % Calculate minimum pair
    minChord = Inf;
    theseChords = ThisConfig.chordList;
    for itrC = 1:length(theseChords)
        thisMin = min(theseChords{itr});
        if (thisMin < minChord)
            minChord = thisMin;
        end
    end
    minPair(itr) = thisMin;
    
    % Calculate mean pair
    meanPair(itr) = mean( [theseChords{:}] );
    
    % Calculate minimum distance to boundary
    theseCoords = ThisConfig.coordVec;
    minDist = Inf;
    for itrD = 1:N
        [~, thisDist] = findDist2Bound(ecc, theseCoords(itrD, :), 0);
        if (thisDist < minDist)
           minDist = thisDist; 
        end
    end
    minBound(itr) = minDist;
    
end


% Remove entries which weren't calculated
dudInd = find(N_vals == 0);
N_vals(dudInd) = [];
distPerTrap(dudInd) = [];
minPair(dudInd) = [];
meanPair(dudInd) = [];
minBound(dudInd) = [];


%  Plot results
% figure(1);
% plot(minPair, 2*minBound, '- .');
% xlabel('Minimum Pair-wise Distance');
% ylabel('Effective Diameter at Border');
% grid on
% 
figure(1);
hold off
plot(N_vals, distPerTrap, '- .k');
hold on
plot(N_vals, meanPair, '- *k', 'MarkerSize', markSize);
plot(N_vals, minPair, '- vk', 'MarkerSize', markSize)
plot(N_vals, 2*minBound, '- ok', 'MarkerSize', markSize)
hold off
legend('Measure of Area per Trap', 'Mean Mutual Distance', 'Minimum Mutual Distance',  '2\cdot Minimum Distance to Boundary');
% legend('2/sqrt(N)', 'mean( |x_n - x_k| )', 'min( |x_n - x_k| )', '2\cdotmin(1 - r)');
xlabel('Number of Traps');

end