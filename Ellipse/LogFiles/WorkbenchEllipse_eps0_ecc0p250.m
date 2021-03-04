clear;

addpath('Auxilliary Code');
addpath('Ellipse');


% User Input
% Domain setup
ecc = 0.25;                  % Eccentricity of ellipse
epsilon = 0;                 % Trap size
% Sweep parameters
N_list = [1:50];     % Trap nums to optimize
diff_thresh = 1e-3;  
optLimit = 1e3;      % Number of optimization iterations
% File name
log_dir = ['Ellipse\LogFiles\Log_eps0-ecc0p', num2str(1000*ecc, 3)];
compData = {[NaN, NaN]}; % Dummy data
compData = readLogFile([log_dir, '\GoodLog21.csv']); % Current best data
compN = getN(compData);
flag_init = 0; % Generate first set of optimums, rather than read from file

% Log file init
log_filename_base = 'log';
comp_log_filename_base = [log_dir, '\', 'GoodLog'];
file_ext = '.csv';

% If the log directory does not exist, create it
if isempty(dir(log_dir))
    mkdir(log_dir);
end

InputStruct = buildInputStruct(ecc, epsilon);

% Optimization parameters
opts = optiset('maxiter', 1e9, 'maxfeval', 1e9, ...
    'tolrfun', 1e-5, 'tolafun', 0, ...
    'display', 'iter');
opts.solverOpts = pswarmset('swarm_size', 100);


%% Run Optimization

a = InputStruct.a;
b = InputStruct.b;
A = InputStruct.A;
mu = InputStruct.mu;

% Initialization
nN = length(N_list);
status_list = zeros(1, nN);

local_search = 'local';
global_search = 'global';


optCount = 0; % Number of optimization iterations
prev_search = global_search;

N_list_old = [];
nN_old = length(N_list_old);
nN = length(N_list);
while ( ~isempty(N_list) && (optCount < optLimit) )
    optCount = optCount+1;
    
    if (strcmpi(prev_search, local_search))
        disp('Global search');
        prev_search = global_search;
    elseif (strcmpi(prev_search, global_search))
        disp('Local search');
        prev_search = local_search;
    end
    
    [ fileID, log_filename, log_filepath, fileInd ] = createLogfile(log_dir, log_filename_base, file_ext);
    % Optimize Configurations
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
                % Distribute traps in elliptical coordinates to avoid
                % clustering
                x0 = buildInitialCoords(A, mu, N);
            else
                thisInd = find( compN == N );
                x0 = compData{thisInd}(2:nvars+1); % Use best known coords
            end
            % Min/Max coordinate bounds
            [ ub, lb ] = buildBounds(A, N);
            
            % Run search algorithm
            if (strcmpi(prev_search, local_search))
                [xOpt, fval, exitflag, info] = opti_pswarm(meritFunc, lb, ub, x0, [], [], opts);
            elseif (strcmpi(prev_search, global_search))
                [xOpt, fval, exitflag, info] = fminsearch(meritFunc, x0);
            end
            status_list(itr) = exitflag;
            
            % Convert to row vector
            [xR, ~] = size(xOpt);
            if (xR ~= 1)
                xOpt = xOpt';
            end
            % Convert to cartesian coords
            [xC, yC] = calcElliptical2Cartesian(xOpt(1:N), mu, xOpt(N+1:end));
            x = [xC, yC];
            
            % Update debugging info
            CompTime = toc(CompHandle);
                       
            % Update Log
            fprintf(fileID, '%u,', N);
            fprintf(fileID, '%0.15f,', x);
            fprintf(fileID, '%0.15f,', fval);
            fprintf(fileID, '%0.15f\n', CompTime);
            
        end
        
    catch ME
        % Close file if an error occurs
        fclose(fileID);
        rethrow(ME);
    end
    fclose(fileID);
    
    % Initialization stage over
    if (flag_init)
       % Check initial set for bad optimums. If any are found, try to
       % improve them
       goodData = readLogFile(log_filepath);
       compData = goodData;
       N_list = identifyBadConfigurations(compData, goodData, diff_thresh);
       flag_init = 0; 
    else 
       goodData = readLogFile([comp_log_filename_base, num2str(fileInd-1), file_ext]); % Load most recent archive
       compData = readLogFile(log_filepath);
       N_list = identifyBadConfigurations(compData, goodData, diff_thresh);
    end
    compN = getN(compData); % Used to map N_list to compData
    nN = length(N_list);
    
    % Write current optimums to log file
    [ writeData ] = updateDataSet(compData, goodData);
    writeLogFile(writeData, [comp_log_filename_base, num2str(fileInd), file_ext]);
    
end
