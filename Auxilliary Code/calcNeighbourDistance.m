function [ chordList ] = calcNeighbourDistance(thisN, coordVec, triD)
% Author: Jason Gilbert
% Data: March 19, 2020
% Version: V00
% Last Updated: N/A
% 
% Input:
%   thisN - number of traps
%   coordVec - Mx2 array of trap coords [x, y]
%   triD - array of indices defining Delaunay triangulation of traps
%          - see documentation for delaunay()
%
% Output:
%   chordList - a cell array of distances to nearest neighbours, calculated
%               using Delaunay triangulation, formatted such that
%               coordVec(ind, :) -> chordList{ind}


% Calculate pair-wise distances
[nTri, mTri] = size(triD);
nChords = 1/2*mTri*(mTri - 1); % Number of chords in each triangle
chords = zeros(nTri, nChords); % Preallocate array to hold chords in each triangle


% Preallocate chord list
chordList = cell(1, thisN);
for itrN = 1:thisN
    nInst = 2*nnz(triD == itrN); % Count number of triangles trap appears in
    chordList{itrN} = zeros(1, nInst); % Add one for image trap
end


% Loop through triangles
for itrN = 1:nTri
    thisTri = triD(itrN, :);
    theseTraps = coordVec(thisTri, :);
    ind = 0;
    
    % Loop through triangle vertices (pairs of traps)
    for itrM = 1:mTri-1
        trapM = theseTraps(itrM, :);
        for itrK = itrM+1:mTri
            ind = ind+1;
            trapK = theseTraps(itrK, :);
            
            thisR = sqrt( sum((trapM - trapK).^2) ); % Chord length
            chords(itrN, ind) = thisR;
            
            % Add chord length to list (fill first empty entry)
            newInd = find( chordList{thisTri(itrM)} == 0, 1 );
            chordList{thisTri(itrM)}(newInd) = thisR;
            % Chord m -> k is the same as k -> m
            newInd = find( chordList{thisTri(itrK)} == 0, 1 );
            chordList{thisTri(itrK)}(newInd) = thisR;
        end
    end
end


% Remove duplicate chords, and sort by chord length
% Each trap may appear in multiple triangles
for itr = 1:thisN
    chordList{itr} = unique(chordList{itr});
end

end