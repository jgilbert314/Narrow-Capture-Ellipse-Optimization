function [  ] = plotDense(hist2D, bins, filter, thresh)
% Author: Jason Gilbert 
% Date: April 25, 2020
% Last Updated: N/A
% Version: V00
% 
% Summary:
%   This function produces a 2D density plot of the radial positions, and 
%   the distances to the nearest neighbours, of the traps in a 
%   configuration. The plot is produced by applying a user-defined filter
%   to a 2D histogram via 2D convolution.
% 
% Input:
%   hist2D - a MxM histogram where the rows represent radial bins, and 
%            columns represent pair-wise bins
%   bins   - a vector defining the bin edges associated with the histogram
%   filter - a vector defining the filter to be applied to the histogram
%            via 2D convolution
%   thresh - colour threshold. Above this threshold, all values are
%            assigned the same colour. 
%            Max: 1
% 
% Output:
%   Overwrites current figure


lenWin = length(filter);

% Calculate trap density vs radii and pair-wise distance
dense2D = conv2(filter, filter, hist2D);
dense2D = dense2D/max(max(dense2D));
dense2D = dense2D(lenWin/2:end-lenWin/2, lenWin/2:end-lenWin/2);
% Calculate normalized projection onto radial axis
denseR = sum(dense2D, 2);
denseR = denseR/max(denseR);
% Calculate normalized projection onto distance axis
denseD = sum(dense2D, 1);
denseD = denseD/max(denseD);


% Plot 2D density
subplot(2, 2, 3);
pcolor(bins, bins, dense2D');
caxis([0, thresh]);
shading flat;
xlabel('Radii');
ylabel('Distances');


% Plot radial density
subplot(2, 2, 1);
plot( bins, denseR, '-');
ylabel({'Radial Density', '(Normalized)'});
grid on;

% Plot pair-wise density
subplot(2, 2, 4);
plot(denseD, bins);
pltTitle = title({'Pair-wise Density', '(Normalized)'});
set(pltTitle, 'FontWeight', 'normal'); % Un-bold title
grid on;


end