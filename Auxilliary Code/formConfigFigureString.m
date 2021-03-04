function [ output_string ] = formConfigFigureString(file_prefix, file_suffix, header)
% Inputs:
%   file_prefix - string, filepath and file prefix common to each file
%   file_suffix - string, portion of filename unique to this figure
%   header - string, comment identifying figure


output_string = { ... 
['% ' header], ...
'\begin{figure}[H]', ...
'\begin{subfigure}[]{', ...
['\includegraphics[width=0.5\textwidth]{Figures/Config/TrapConfig_ecc0p0_N', file_suffix, '}'], ...
'}', ...
'\end{subfigure}', ...
'\begin{subfigure}[]{', ...
['\includegraphics[width=0.5\textwidth]{Figures/Config/TrapConfig_ecc0p125_N', file_suffix, '}'], ...
'}', ...
'\end{subfigure}', ...
'\begin{subfigure}[]{', ...
['\includegraphics[width=0.5\textwidth]{Figures/Config/TrapConfig_ecc0p250_N', file_suffix, '}'], ...
'}', ...
'\end{subfigure}', ...
'\begin{subfigure}[]{', ...
['\includegraphics[width=0.5\textwidth]{Figures/Config/TrapConfig_ecc0p500_N', file_suffix, '}'], ...
'}', ...
'\end{subfigure}', ...
'\caption{***}', ...
'\end{figure}', ...
};

output_string = sprintf('%s\n', output_string{:});

end