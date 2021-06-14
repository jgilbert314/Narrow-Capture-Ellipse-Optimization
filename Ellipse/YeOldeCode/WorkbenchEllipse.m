%% ATTENTION: THIS SCRIPT IS OBSOLETE. Updated versions use optimizationOverseer()
clear;

addpath('Auxilliary Code');
addpath('Ellipse');


% User Input
% Domain setup
ecc = 0;                  % Eccentricity of ellipse
epsilon = 0.05;                 % Trap size
% Sweep parameters
N_list = [1:50];     % Trap nums to optimize
diff_thresh = 1e-3;  
optLimit = 1e3;      % Number of optimization iterations
% File name
log_dir = ['Ellipse\LogFiles\Log_eps0p500-ecc0p', num2str(1000*ecc, 3)];
compData = {[NaN, NaN]}; % Dummy data
compData = readLogFile([log_dir, '\GoodLog235.csv']); % Current best data
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
opts = optiset('maxiter', 1e6, 'maxfeval', 1e6, ...
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
circ_flag = 0;
if (a == b) % Circular domain
   circ_flag = 1; 
end

local_search = 'local';
global_search = 'global';


optCount = 0; % Number of optimization iterations
prev_search = global_search;


N_list_old = [];
nN_old = length(N_list_old);
nN = length(N_list);
while ( ~isempty(N_list) && (optCount < optLimit) )
    optCount = optCount+1;
    
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
            meritFunc = @(x) meritFuncGeneral(x, InputStruct);
            
            % Update init coords and bounds
            if (flag_init)
                % Distribute traps in elliptical coordinates to avoid
                % clustering
                x0 = buildInitialCoords(A, mu, N);
            else
                thisInd = find( compN == N );
                x0 = compData{thisInd}(2:nvars+1); % Use best known coords
                if (circ_flag)
                   [r, theta] = calcCartesian2Polar(x0(1:N), x0(N+1:end));
                   x0 = [r, theta];
                else % NOTE: untested
                   [A_p, nu_p] = calcCartesian2Elliptical(x0(1:N), x0(N+1:end), mu);
                   x0 = [A_p, nu_p];
                end
            end
            % Min/Max coordinate bounds
            [ ub, lb ] = buildBounds(A, N);
            
            % Run search algorithm
            if (strcmpi(prev_search, local_search))
                [xOpt, fval, exitflag, info] = opti_pswarm(meritFunc, lb, ub, x0, [], [], opts);
                prev_search = global_search;
            elseif (strcmpi(prev_search, global_search))
                [xOpt, fval, exitflag, info] = fminsearch(meritFunc, x0);
                prev_search = local_search;
            end
            status_list(itr) = exitflag;
            
            % Convert to row vector
            [xR, ~] = size(xOpt);
            if (xR ~= 1)
                xOpt = xOpt';
            end
            % Convert to cartesian coords
            if (circ_flag)
                [xC, yC] = calcPolar2Cartesian(xOpt(1:N), xOpt(N+1:end));
            else
                [xC, yC] = calcElliptical2Cartesian(xOpt(1:N), mu, xOpt(N+1:end));
            end
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
