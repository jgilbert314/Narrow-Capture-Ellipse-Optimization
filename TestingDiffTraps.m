%%
% clear;

addpath('Auxilliary Code');
addpath('Ellipse');


% User Input
% Domain setup
ecc = 0.125;                    % Eccentricity of ellipse
epsilon = 0.05;             % Trap size
N_p = 2;                    % Number of traps differing in size
k = 1/2;                    % Size of one trap subset relative to epsilon_0
InputStruct = buildInputStruct(ecc, epsilon); % Initialize arguments to pass to optimizer
InputStruct.epsilon_0 = epsilon; % TESTING
InputStruct.N_p = N_p;
InputStruct.k = k;

% Sweep parameters
InputStruct.N_list = [N_0];     % Trap nums to optimize
InputStruct.diff_thresh = 1e-6;  % Improvements smaller than this value are considered converged
InputStruct.optLimit = 100;      % Maximum number of optimization iterations
InputStruct.init_search = 'local'; % Search algorithm to run on first pass ('global' or 'local')
InputStruct.flag_init = 0;       % Generate first set of optimums, rather than read from file
% File name
InputStruct.log_dir = ['Ellipse\LogFiles\Log_eps-diff', num2str(1000*epsilon, 3), '0-ecc0p', num2str(1000*ecc, 3)];
InputStruct.archive_base = 'GoodLog';
InputStruct.ref_filename = [InputStruct.log_dir, '\GoodLog0.csv'];
% Log file init
InputStruct.log_filename_base = 'log';
InputStruct.comp_log_filename_base = [InputStruct.log_dir, '\', 'GoodLog'];
InputStruct.file_ext = '.csv';


% Optimization algorithm parameters
OptStruct = optiset('maxiter', 1e6, 'maxfeval', 1e6, ...
    'tolrfun', 1e-6, 'tolafun', 0, ...
    'display', 'iter');
OptStruct.solverOpts = pswarmset('swarm_size', 100);
InputStruct.OptStruct = OptStruct;

% Local optimization parameters (as used by fminsearch)
InputStruct.OptStructLocal = optimset('MaxFunEvals', OptStruct.maxfeval, 'MaxIter', OptStruct.maxiter, ...
   'TolFun', OptStruct.tolrfun, 'Display', 'iter', 'TolX', 1e-6);


InputStruct.eps_inds = eps_inds; % TESTING

%% Run optimization

optimizationOverseer(InputStruct);

