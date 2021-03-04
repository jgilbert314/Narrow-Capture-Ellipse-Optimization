clear;

addpath('Auxilliary Code');
addpath('Ellipse');


% User Input
N_list = [1:100];
diff_thresh = 1e-3;
log_dir = 'LogFiles';
refData = {[2, NaN, NaN, NaN, NaN, 0, NaN]}; % Dummy input
flag_init = 1; % Generate first set of optimums, rather than read from file

% Log file init
log_filename_base = 'log';
running_log_filename_base = [log_dir, '\', 'GoodLog'];
file_ext = '.csv';


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


% Micellaneous terms
D = 1;
nu = -1/log(epsilon);

% Input Structure setup
InputStruct = struct();

% Miscellaneous parameters
InputStruct.epsilon = epsilon;
InputStruct.D = D;

% Domain specific parameters
InputStruct.a = a; InputStruct.b = b;
InputStruct.f = f; 
InputStruct.beta = beta; InputStruct.xi_b = xi_b;
InputStruct.omega = omega;


% Optimization parameters
opts = optiset('maxiter', 1e9, 'maxfeval', 1e9, ...
    'tolrfun', 1e-5, 'tolafun', 0, ...
    'display', 'iter');
opts.solverOpts = pswarmset('swarm_size', 100);


%% Run Optimization

% Initialization
nN = length(N_list);
opt_param = cell(1, nN);
p = zeros(1, nN);
status_list = zeros(1, nN);
CompTime = zeros(1, nN);

% If the log directory does not exist, create it
if isempty(dir(log_dir))
    mkdir(log_dir);
end


optCount = 0; % Number of optimization iterations

N_list_old = [];
nN_old = length(N_list_old);
nN = length(N_list);
while ( ~isempty(N_list) && (optCount < 1e3) )
    optCount = optCount+1;
    
    % Create new log file
    fileInd = 0;
    log_filename = [log_filename_base, num2str(fileInd), file_ext];
    while (checkFileExists(log_dir, log_filename))
        fileInd = fileInd+1;
        log_filename = [log_filename_base, num2str(fileInd), file_ext];
    end
    log_filepath = [log_dir, '\', log_filename];
    
    % Optimize Configurations
    fileID = fopen(log_filepath, 'w'); % Create log file
    try
        for itr = 1:nN
            CompHandle = tic;
            
            N = N_list(itr);
%             opts.solverOpts = pswarmset('swarm_size', N);
            nvars = 2*N;
            
            % Update merit function
            InputStruct.N = N; 
            meritFunc = @(x) meritFuncEllipse(x, InputStruct);
            
            % Update init coords and bounds
            if (flag_init)
                x0 = [ (rand(1, N) - 0.5)*a, (rand(1, N) - 0.5)*b]'; % Initial config
            else
                thisInd = find( compN == N );
%                 x0 = compData{thisInd}(2:nvars+1); % Use current optimal coords
                x0 = (1/sqrt(2))*rand(1, nvars);
            end
            % Max coordinate bounds
            B = [a*ones(1, N), b*ones(1, N)]';
            lb = -B; ub = B; 
            
            [x, fval, exitflag, info] = opti_pswarm(meritFunc, lb, ub, x0, [], [], opts);
            status_list(itr) = exitflag;
            
            opt_param{itr} = x;
            p(itr) = fval;
            thisTime = toc(CompHandle);
            CompTime(itr) = thisTime;
            
            outputLine = [N, x', fval, thisTime]; % x is output as 2Nx1
            
            % Update Log
            fprintf(fileID, '%u,', N);
            fprintf(fileID, '%0.15f,', x);
            fprintf(fileID, '%0.15f,', fval);
            fprintf(fileID, '%0.15f\n', CompTime(itr));
            
        end
        
    catch ME
        fclose(fileID);
        rethrow(ME);
    end
    fclose(fileID);
    
    % Initialization stage over
    if (flag_init)
       goodData = readLogFile(log_filepath);
       flag_init = 0; 
    else 
       goodData = readLogFile([running_log_filename_base, num2str(fileInd-1), file_ext]);
    end
    
    % Compare new optimums to reference. If improvement is necessary,
    % the loop will repeat
    compData = readLogFile(log_filepath);
    N_list = identifyBadConfigurations(goodData, refData, diff_thresh);
    compN = getN(compData); % Used to map N_list to compData
    nN = length(N_list);
    
    % Write current optimums to log file
    [ writeData ] = updateDataSet(compData, goodData);
    writeLogFile(writeData, [running_log_filename_base, num2str(fileInd), file_ext]);
    
end
