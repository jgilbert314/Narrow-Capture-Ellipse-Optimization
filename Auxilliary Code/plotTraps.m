function [] = plotTraps(trapCoords, varargin)
% Author: Jason Gilbert
% Date: March 22, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function plots the boundary of the domain, and trap coordinates 
%   in one subplot, and the radial coordinates of the traps in another.
%   Optionally, a Delaunay triangulation of the trap coordinates may be
%   added to the first subplot, if specified by name-value pair.
% 
% Input:
%   trapCoords - an Nx2 matrix or trap coordinates, with rows in [x, y]
%                format
%   varagin - optional arguments. Currently supported:
%       tri - if one of the optional arguments is the string 'tri'
%             (case insensitive), and the next is an Mx3 matrix of 
%             indices produced by Delaunay triangulations (see delaunay() 
%             documentation), the triangulation will be plotted.


% Optional argument handling
flagTri = 0; % Flag that triangulation should be plotted
if ( ~isempty(varargin) ) % Check if optional arguements given
    numVarargs = length(varargin);
    if (  mod(numVarargs, 2) ) % Check if all name-value pairs are complete
       error('Incomplete name-value pair given,'); 
    end
    
    % Check if valid argument has been given
    argFlag = 0;
    varargList = {'tri'};
    for itr = 1:numVarargs
        if ( any(strcmpi(varargList{itr}, varargin)) )
            argFlag = 1;
            break;
        end
    end
    if (~argFlag)
        error('No valid optional arguments');
    end
    
    % Check if triangulation should be plotted. If so, extract 
    % triangulation data
    if ( any(strcmpi('tri', varargin)) )
        flagTri = 1;
        triInd = find(strcmpi('tri', varargin));
        trapTri = varargin{ triInd+1 };
    end
end



% Extract trap coordinates for plotting
tX = trapCoords(:, 1);
tY = trapCoords(:, 2);
r = sqrt( sum( trapCoords.^2, 2 ) );
r = sort(r);
thisN = length(r);

% Boundary coords
bX = linspace(-1, 1, 1e3);
bY = sqrt( 1 - bX.^2 );


figure();
subplot(2, 1, 1);
hold off
plot(tX, tY, '.k', 'MarkerSize', 15); % Plot trap coords
hold on
if (flagTri) % Plot trap triangulation
    triplot(trapTri, tX, tY, 'k'); % Plot triangulation
end
plot(bX,  bY, 'k');  % Plot upper boundary
plot(bX, -bY, 'k');  % Plot lower boundary
hold off
title({'Trap Positions' ; ['N = ', num2str(thisN)] });
grid on
axis equal

% Plot radial coordinates of traps
subplot(2, 1, 2);
hold off
plot(r, '- .k');
hold on
ylabel('Radial Coordinates');
ylim([0, 1]);
grid on


end