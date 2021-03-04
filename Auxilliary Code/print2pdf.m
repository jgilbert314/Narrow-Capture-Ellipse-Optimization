function [  ] = print2pdf(F, filename)
% Author: Jason Gilbert
% Date: August 12, 2020
% Version: V00
% Last Updated: N/A
% 
% Summary:
%   This function saves a figure as a pdf, sized to accomodate the figure
% 
% Input:
%   F - handle of figure to be saved as PDF
%   filename - name of file to be produced

set(F, 'Units', 'Inches');
pos = get(F, 'Position');
set(F, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
print([filename, '.pdf'], '-dpdf', '-fillpage');

end