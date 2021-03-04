%% Initialization
clear;

% Plotting setup 
save_flag = 0; % Flags that all data should be visualized and saved
SaveDir = ['VisFiles']; % Relative path to plot save folder
% Domain setup
filename_list = {'Ellipse\LogFiles\Log_eps0p500-ecc0p500\log1.csv'};
ecc_list = [0.500];

% filename = 'Ellipse\LogFiles\Log_eps0p500-ecc0p500\GoodLog3.csv';



%% Compare Image Distance to Max Radial coordinate
for itr = 1:length(filename_list)
    filename = filename_list{itr};
    
    DataSet = buildDataSet(filename);
    plotTrapDist(DataSet, ecc_list(itr));
end


%%
[a, b] = calcAxis(ecc_list(1));

ind = 25;
ThisData = DataSet.ConfigData(ind);
N = ThisData.thisN;
trapCoords = ThisData.coordVec;
plotTrapsEllipse(trapCoords, a, b)