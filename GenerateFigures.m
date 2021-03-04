%% Initialization
clear;
addpath('Results');
addpath('Ellipse');
addpath('Auxilliary Code');
% Specify file path delimiter based on OS
if (ispc) % Windows
   path_delim = '\'; 
elseif (isunix) % Linux
   path_delim = '/';
elseif (ismac) % Mac (untested)
   path_delim = ':';
end

filename_list = {'Log_eps0p500-ecc0p000.csv', 'Log_eps0p500-ecc0p125.csv', 'Log_eps0p500-ecc0p250.csv', 'Log_eps0p500-ecc0p500.csv'};
ecc_list = {0, 0.125, 0.250, 0.500};
N_min = 3; % Minimum number required for triangulation is 3
N_max = 50;


% Global Plot Parameters
FigureType = struct('png', 'png', 'pdf', 'pdf');
figure_type = 'pdf';
grid_flag = 0; % Flags that background grid should be applied
title_flag = 1; % Flags that each plot should be titled
font_size = 16; % Specifies text font


save_flag = 1; % Flag indicating plots should be saved
save_dir = formFilepath( {'Paper'}, path_delim );


num_files = length(filename_list);

%% Plot Trap Configurations

N_list = [5, 10, 25, 40]; % Configurations to plot

distMarkSize = 15; % Default marker size for trap coordinates
radMarkSize = 12;  % Default marker size for radial coordinates

if ( strcmp(figure_type, FigureType.png) )
    sub_dir = '';
elseif ( strcmp(figure_type, FigureType.pdf) )
    sub_dir = '';
else
    error('Unsupported image export format.');
end
this_dir = formFilepath( {save_dir, sub_dir}, path_delim );
mkdir(this_dir);
close('all');

for itr = 1:num_files
    
    DataSet = buildDataSet( filename_list{itr} ); % Initialize configuration data
    ecc = ecc_list{itr};
    [a, b] = calcAxis( ecc );           % Major and minor axis of domain
    
    % Plot trap configurations
    for itrN = 1:length(N_list)
        thisN = N_list(itrN);
        thisCoord = DataSet.ConfigData(thisN).coordVec;
        thisTri = DataSet.ConfigData(thisN).triD;
        plotTrapsEllipse(thisCoord, a, b, 'tri', thisTri, ... 
            'distMarkSize', distMarkSize, 'radMarkSize', radMarkSize);
        
        F = gcf;
        % Apply grid to all figures
        if (grid_flag)
           for itrFig = 1:length(F.Children)
                grid(F.Children(itrFig), 'on'); 
           end
        end
        % Use same scale for all plots
        F.Children(2).XLim = [-3.2, 3.2];
        
        if (save_flag)
            this_filepath = [this_dir, path_delim, 'TrapConfig_ecc0p', num2str(ecc*1000, 3), '_N', num2str(thisN)];
            if (strcmp(figure_type, FigureType.png))
                print( this_filepath, '-r300', '-dpng' );
            elseif (strcmp(figure_type, FigureType.pdf))
                F = gcf;
                print2pdf(F, this_filepath);
            end
        end
        
    end
    
end

%% Plot Merit Function

marker_size = 12.5;

sub_dir = '';
this_dir = formFilepath( {save_dir, sub_dir}, path_delim );
mkdir(this_dir);
close('all');

for itr = 1:num_files
    
    dataArr = readLogFile( filename_list{itr} ); % Initialize configuration data
    dataArr = dataArr(1:N_max);
    ecc = ecc_list{itr};
    
    p = getIntEng(dataArr);
    theseN = getN(dataArr);
    plot(theseN, p, '- . k', 'MarkerSize', marker_size);
    
    
    xlabel('Number of Traps');
    if (title_flag)
        title( ['Eccentricity: ', num2str(ecc)] );
    end
    ylim([-0.25, 0.20])
    xlim([0, N_max])
    if (grid_flag)
        grid on;
    end
    
    if (save_flag)
        this_filepath = [this_dir, path_delim, 'MeritFunc_ecc0p', num2str(ecc*1000, 3)];
        if (strcmp(figure_type, FigureType.png))
            print( this_filepath, '-r300', '-dpng' );
        elseif (strcmp(figure_type, FigureType.pdf))
            F = gcf;
            print2pdf(F, this_filepath);
        end
    end
end

%% Plot Distance Comparison

legend_fontsize = 12;
marker_size = 5;

sub_dir = '';
this_dir = formFilepath( {save_dir, sub_dir}, path_delim );
mkdir(this_dir);
close('all');

for itr = 1:length(filename_list)
    filename = filename_list{itr};
    ecc = ecc_list{itr};
    
    DataSet = buildDataSet(filename);
    DataSet.ConfigData = DataSet.ConfigData(1:N_max);
    DataSet.structIndex = DataSet.structIndex(1:N_max);
    plotTrapDist(DataSet, ecc, marker_size);
    if (title_flag)
        title(['Eccentricity: ', num2str(ecc)]);
    end
    if (grid_flag)
       grid on; 
    end
    fig = gca;
    fig.Legend.FontSize = legend_fontsize;
    
    if (save_flag)
        this_filepath = [this_dir, path_delim, 'TrapDist_ecc0p', num2str(ecc*1000, 3)];
        if (strcmp(figure_type, FigureType.png))
            print( this_filepath, '-r300', '-dpng' );
        elseif (strcmp(figure_type, FigureType.pdf))
            F = gcf;
            print2pdf(F, this_filepath);
        end
    end
    
end

%% Generate Table Strings

% Tabulate values
col_num = 5;
dataTable = zeros(N_max+1, col_num);
dataTable(1, 1) = 0;

for itr = 1:num_files
    
    dataArr = readLogFile( filename_list{itr} ); % Initialize configuration data
    dataArr = dataArr(1:N_max);
    ecc = ecc_list{itr};
    
    p = getIntEng(dataArr);
    theseN = getN(dataArr);
    
    % Add values to table
    dataTable(2:end, itr+1) = p;
    dataTable(1, itr+1) = ecc;
    
end

% Add N values to table
dataTable(2:end, 1) = theseN;

% Write table body to file
bodyString = [];
[num_rows, num_cols] = size(dataTable);
delim = ' & ';          % Latex table element seperator
line_delim = ' \\\\\\\\\n'; % Latex table line seperator (formated for sprintf)
format_spec = [ '%u', repmat([delim, '%0.4f'], [1, num_cols-1]), line_delim ];

for itrR = 1:num_rows
    bodyString = [ bodyString, sprintf(format_spec, dataTable(itrR, :)) ];
end

clc;
fprintf(bodyString);

%% Generate Comparison Figure/Table

marker_size = 12.5;

data_path = 'C:\Projects\Narrow Escape\Unit Disk\Archive\WardData\ref_p_cell';
load(data_path); % Loads variable named ref_p


% Load old data
refN = getN(ref_p);
refP = getIntEng(ref_p);


% Load new data
dataArr = readLogFile( filename_list{1} ); % Load new values
dataArr = dataArr(1:N_max);

newP = getIntEng(dataArr);
newN = getN(dataArr);
newP = newP(refN);
newN = newN(refN);


% Calculate AMFPT and comparison
ecc = ecc_list{1};
D = 1;
eps = 0.05;
omega = pi;
nu = -1/log(eps);
refU = omega./(2*pi*D*nu.*refN) .* (1 + 2*pi*nu./refN.*refP);
newU = omega./(2*pi*D*nu.*newN) + 2*pi./newN.*newP;

relDiff = (newU - refU)./refU;


plot(refN, relDiff, '- . k', 'MarkerSize', marker_size);
xlabel('Number of Traps');
ylabel('Relative Difference between Average MFPT');
if (grid_flag)
    grid on;
end

sub_dir = '';
this_dir = formFilepath( {save_dir, sub_dir}, path_delim );
mkdir(this_dir);
if (save_flag)
    this_filepath = [this_dir, path_delim, 'RelativeDifference'];
    if (strcmp(figure_type, FigureType.png))
        print( this_filepath, '-r300', '-dpng' );
    elseif (strcmp(figure_type, FigureType.pdf))
        F = gcf;
        print2pdf(F, this_filepath);
    end
end


% Write table body to file
dataTable = [ refN' , refU', newU' ];

bodyString = [];
[num_rows, num_cols] = size(dataTable);
delim = ' & ';          % Latex table element seperator
line_delim = ' \\\\\\\\\n'; % Latex table line seperator (formated for sprintf)
format_spec = [ '%u', repmat([delim, '%0.5f'], [1, num_cols-1]), line_delim ];

for itrR = 1:num_rows
    bodyString = [ bodyString, sprintf(format_spec, dataTable(itrR, :)) ];
end

clc;
fprintf(bodyString);