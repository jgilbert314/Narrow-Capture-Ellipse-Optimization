%% Initialization
clear;

addpath('Auxilliary Code');

save_flag = 0; % Flags that all data should be visualized and saved
SaveDir = ['VisFiles']; % Relative path to plot save folder

filename = 'Ellipse\LogFiles\log4.csv';
dataArr = readLogFile(filename);

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

%% Save Plots of Trap Coords for all Configurations

save_flag = 0;

% Domain setup
ecc = 0.3;                % Eccentricity of ellipse
epsilon = 0;              % Trap size

% Domain dependent terms
[a, b] = calcAxis(ecc); % Axis of ellipse
omega = pi*a*b;
% Equation (4.3b)
f = sqrt(a^2 - b^2);
beta = (a - b)/(a + b);
xi_b = -log(beta)/2;
ecc = 0.3;                % Eccentricity of ellipse
epsilon = 0;              % Trap size

% Domain dependent terms
[a, b] = calcAxis(ecc); % Axis of ellipse

res  = 1e3;
x = linspace(-a, a, res);
y = b*sqrt(1 - (x/a).^2);



for itr = 3:length(DataSet.ConfigData)
    thisCoord = DataSet.ConfigData(itr).coordVec;
    thisTri = DataSet.ConfigData(itr).triD;
    close('all');
    plotTrapsEllipse(thisCoord, a, b, 'tri', thisTri);
    subplot(2, 1, 1);
    hold on
    plot(x, y, 'k');
    plot(x, -y, 'k');
    hold off
    axis equal;
    if (save_flag)
        print( [SaveDir, '\N', num2str(itr)], '-r300', '-dpng' );
    end
end

%% Compare Image Distance to Max Radial coordinate

% load('DataSet');

numSets = length(DataSet.ConfigData);
minPair = zeros(1, numSets);
minBound = zeros(1, numSets);

for itr = 3:numSets
    ThisConfig = DataSet.ConfigData(itr);
    
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
    
    % Calculate minimum boundary
    theseCoords = ThisConfig.coordVec;
    tX = theseCoords(:, 1);
    tY = theseCoords(:, 2);
    r = sqrt( (tX/a).^2 + (tY/b).^2 );
    minBound(itr) = min( 1 - r );
    
end

maxR = 1 - minBound;
imgDist = 1./maxR - maxR;

% Remove entries for N=1, N=2
minPair = minPair(3:end);
imgDist = imgDist(3:end);

figure(1);
plot(minPair, imgDist, '- *');
xlabel('Minimum Pair-wise Distance');
ylabel('Distance to Image Trap');
grid on

figure(2);
hold off
plot(minPair, '- *k')
hold on
plot(2*minBound, '- .')
hold off
legend('min( |x_n - x_k| )', '2\cdotmin(1 - r)');
xlabel('Number of Traps');
grid on

%% Demo Density Plot

demo_flag = 1;
res = 1e3;
thresh = 0.3;

% Init Demo Visualization
ind = 20; % Index of demo data
ThisConfig = DataSet.ConfigData(ind);
thisTri = ThisConfig.triD;
thisCoord = ThisConfig.coordVec;
thisChordList = ThisConfig.chordList;
thisN = ThisConfig.thisN;

tX = thisCoord(:, 1);
tY = thisCoord(:, 2);
r = sqrt( sum( thisCoord.^2, 2 ) );
[r, rInd] = sort(r);                  % Sort radii from least to greatest
thisChordList = thisChordList(rInd);  % Sort chords by the radii of central trap


% Build base histogram
edges = linspace(0, 1, res+1);
[hist, ~, histInd] = histcounts(r, edges);
[ hist2D, bins ] = distHist(res, r, thisChordList);

% Build filter
resT = res;
t = linspace(-1, 1, resT);
sigma = 0.005;
window = exp( -t.^2/2/sigma^2 ) / sqrt(2*pi) * sqrt(sigma); % Gaussian Filter


if (demo_flag)
    
    % Plot Raw Hist
    figure(1)
    subplot(2, 1, 1);
    bar(bins, hist);
    title('Radial Histogram');
    
    subplot(2, 1, 2);
    pcolor(bins, bins, (hist2D'));
    shading flat;
    caxis([0, thresh]);
    title('Raw 2D Histogram');
    
    % Plot Filtered Density
    figure(2);
    plotDense(hist2D, bins, window, 0.2);
    colormap(jet);
    
end


%% Calculate and plot densities

SaveDir = 'VisFiles';

res = 1e3; % Resolution of density plot

% Build filter
resT = res;
t = linspace(-1, 1, resT);
sigma = 0.005;
window = exp( -t.^2/2/sigma^2 ) / sqrt(2*pi) * sqrt(sigma); % Gaussian Filter

for itr = 3:length(DataSet.ConfigData)
    
    ThisConfig = DataSet.ConfigData(itr);
    thisTri = ThisConfig.triD;
    thisCoord = ThisConfig.coordVec;
    thisChordList = ThisConfig.chordList;
    thisN = ThisConfig.thisN;
    
    tX = thisCoord(:, 1);
    tY = thisCoord(:, 2);
    r = sqrt( sum( thisCoord.^2, 2 ) );
    [r, rInd] = sort(r);                  % Sort radii from least to greatest
    thisChordList = thisChordList(rInd);  % Sort chords by the radii of central trap
    
    
    % Calculate Base Histogram
    [ hist2D, bins ] = distHist(res, r, thisChordList);
    % Generate Plot
    close('all');
    colormap(jet)
    plotDense(hist2D, bins, window, 0.1);
    subplot(2, 2, 1);
    title(['N = ', num2str(thisN)]);
    if (save_flag)
        print( [SaveDir, '\Dense_N', num2str(itr)], '-r300', '-dpng' );
    else
        drawnow();
        pause(0.1);
    end
    
    
end


