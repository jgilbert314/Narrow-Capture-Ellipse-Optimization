function [ DataSet ] = buildDataSet(log_filename)
% *** COMMENT

dataArr = readLogFile(log_filename);

N = getN(dataArr);
numN = length(N);
ConfigData(numN) = struct('thisN', [], 'coordVec', [], 'chordList', [], 'triD', []);
structIndex = [ N ; zeros(1, numN) ]; % [ N ; structure index]
DataSet = struct('ConfigData', ConfigData, 'structIndex', structIndex);

for itr = 1:numN
    thisN = N(itr);
    if (thisN > 2) % Triangulation undefined for N < 2
        DataSet.structIndex(:, itr) = [thisN ; itr];
        ThisConfig = struct('thisN', [], 'coordVec', [], 'chordList', []);
        ThisConfig.thisN = thisN;
        
        
        % Trap Coords
        coords = dataArr{itr}(2:2*thisN+1);
        tX = coords(1:thisN);
        tY = coords(thisN+1:end);
        ThisConfig.coordVec = [ tX', tY' ];
        
        % Triangulate points
        ThisConfig.triD = delaunay(tX, tY);
        ThisConfig.chordList = calcNeighbourDistance(ThisConfig.thisN, ThisConfig.coordVec, ThisConfig.triD);
        
        DataSet.ConfigData(itr) = ThisConfig;
        
    end
end

end