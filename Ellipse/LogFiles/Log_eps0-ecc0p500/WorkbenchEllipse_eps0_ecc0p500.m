%%
clear;

addpath('Auxilliary Code');
addpath('Ellipse');


% User Input
% Domain setup
ecc = 0.500;                  % Eccentricity of ellipse
epsilon = 0;                  % Trap size
InputStruct = buildInputStruct(ecc, epsilon); % Initialize arguments to pass to optimizer


% Sweep parameters
InputStruct.N_list = [1:50];     % Trap nums to optimize
InputStruct.diff_thresh = 1e-6;  % Improvements smaller than this value are considered converged
InputStruct.optLimit = 1e3;      % Maximum number of optimization iterations
InputStruct.init_search = 'local'; % Search algorithm to run on first pass ('global' or 'local')
InputStruct.flag_init = 0;       % Generate first set of optimums, rather than read from file
% File name
InputStruct.log_dir = ['Ellipse\LogFiles\Log_eps0-ecc0p', num2str(1000*ecc, 3)];
InputStruct.archive_base = 'GoodLog';
InputStruct.ref_filename = [InputStruct.log_dir, '\GoodLog23.csv'];
% Log file init
InputStruct.log_filename_base = 'log';
InputStruct.comp_log_filename_base = [InputStruct.log_dir, '\', 'GoodLog'];
InputStruct.file_ext = '.csv';


% Optimization algorithm parameters
OptStruct = optiset('maxiter', 1e9, 'maxfeval', 1e9, ...
    'tolrfun', 1e-5, 'tolafun', 0, ...
    'display', 'iter');
OptStruct.solverOpts = pswarmset('swarm_size', 100);
InputStruct.OptStruct = OptStruct;






%% Run optimization

optimizationOverseer(InputStruct);

