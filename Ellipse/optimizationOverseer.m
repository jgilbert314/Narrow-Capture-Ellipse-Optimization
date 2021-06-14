function [] = optimizationOverseer(InputStruct)

%% Initialization

% Optimization parameters
N_list = InputStruct.N_list;            % Trap nums to optimize
diff_thresh = InputStruct.diff_thresh;  % Improvements smaller than this value are considered converged
optLimit = InputStruct.optLimit;        % Maximum number of optimization iterations
init_search = InputStruct.init_search;  % Search algorithm to run on first pass ('global' or 'local')
flag_init = InputStruct.flag_init;      % Indicates if reference data has been provided for initial guesses
OptStruct = InputStruct.OptStruct;      % Structure of optimization algorithm options

eps_inds = InputStruct.eps_inds;

% Domain parameters
A = InputStruct.A;   % Size in elliptical coordinates
mu = InputStruct.mu; % Eccentricity in elliptical coordinates

% Trap size parameters
epsilon_0 = InputStruct.epsilon_0;
N_p = InputStruct.N_p;
k = InputStruct.k;

% Log file init
log_dir = InputStruct.log_dir; % Directory where log files will be stored
log_filename_base = InputStruct.log_filename_base; % Common log file prefix
ref_filename = InputStruct.ref_filename; % Name of file containing reference data
comp_log_filename_base = InputStruct.comp_log_filename_base; % Common prefix for cumulative log files
file_ext = InputStruct.file_ext; % Common file extension for all output files



% Initialization
nN = length(N_list);
status_list = zeros(1, nN); % Used for debugging optimizer exit status

% Variables which determine the algorithm used
local_search = 'local';
global_search = 'global';
search_type = init_search;

% Check if domain is circular
circ_flag = 0;
if (A == 0)
   circ_flag = 1;
   disp('Domain is circular');
end


% Local optimization parameters (as used by fminsearch)
OptStructLocal = InputStruct.OptStructLocal;


% If the log directory does not exist, create it
if isempty(dir(log_dir))
    mkdir(log_dir);
end

% Check for initialization data
if (~flag_init)
    compData = readLogFile(ref_filename); % Current best data
    compN = getN(compData);
else % Dummy input
    compData = {[NaN, NaN]}; 
    compN = 1:length(compData);
end


%% Optimization
optCount = 0; % Number of optimization iterations
nN = length(N_list);
while ( ~isempty(N_list) && (optCount < optLimit) )
    optCount = optCount+1;
    
    % Debugging messages
    if (strcmpi(search_type, global_search))
        disp('Global Search');
    elseif (strcmpi(search_type, local_search))
        disp('Local Search');
    end
    
    [ fileID, log_filename, log_filepath, fileInd ] = createLogfile(log_dir, log_filename_base, file_ext);
    % Optimize Configurations
    try
        for itr = 1:nN
            CompHandle = tic;
            
            N = N_list(itr);
            nvars = 2*N;
            
            % Update merit function
            InputStruct.N = N;
            InputStruct.epsilon = calcTwoTrapSizes(epsilon_0, N, N_p, k); % TESTING
            InputStruct.epsilon = InputStruct.epsilon(eps_inds); % Shuffle trap positions
            InputStruct.nu = -1./log(InputStruct.epsilon);
            meritFunc = @(x) meritFuncGeneral(x, InputStruct);
            
            % Update init coords and bounds
            if (flag_init)
                % Distribute traps in elliptical coordinates to avoid clustering
                x0 = buildInitialCoords(A, mu, N);
            else
                x0 = compData{compN == N}(2:nvars+1); % Use best known coords
                if (circ_flag) % Convert to polar coordinates
                   [r, theta] = calcCartesian2Polar(x0(1:N), x0(N+1:end));
                   x0 = [r, theta];
                else % Convert to elliptical coords
                    [A_0, nu_0] = calcCartesian2Elliptical(x0(1:N), x0(N+1:end), mu);
                    x0 = [A_0, nu_0];
                end
                
                
            end
            % Min/Max coordinate bounds
            [ ub, lb ] = buildBounds(A, N);
            
            % Run search algorithm
            if (strcmpi(search_type, global_search))
                [xOpt, fval, exitflag, info] = opti_pswarm(meritFunc, lb, ub, x0, [], [], OptStruct);
            elseif (strcmpi(search_type, local_search))
                [xOpt, fval, exitflag, info] = fminsearch(meritFunc, x0, OptStructLocal);
            end
            if (strcmpi(search_type, 'gradient')) % DEBUG - Gradient descent method
                InputStruct.GradStruct.dx = InputStruct.dx*ones(1, nvars);
                [xOpt, exitflag] = gradDescent(meritFunc, x0, InputStruct.GradStruct);
                fval = meritFunc(xOpt);
                info = struct();
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
            
            % Update timing info
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
    
    
    % Update search algorithm
    if (strcmpi(search_type, global_search))
        search_type = local_search;  % Next search will be local
    elseif (strcmpi(search_type, local_search))
        search_type = global_search; % Next search will be global
    end
    
    % TESTING
    indFileID = fopen([comp_log_filename_base, num2str(fileInd), '.ind'], 'w');
    fprintf(indFileID, '%u,', eps_inds); 
    % END TESTING
    
end