function [ x, exitFlag ] = gradDescent(meritFunc, x0, OptStruct)

dx = OptStruct.dx;              % Interval used to calculate gradient
thresh = OptStruct.thresh;      % Threshold of convergence (relative)
iterLim = OptStruct.iterLim;    % Maximum number of iterations

num_params = length(x0); % Number of parameters to optimize

% gamma = zeros(1, num_params); % Vector of parameter increments
gamma = dx;
F0 = meritFunc(x0);           % Function value at initial guess (used to measure convergence)
delF0 = zeros(1, num_params); % Gradient at previous query point

threshFlag = 1; % Flag indicating convergence has been reached
iterFlag = 1;   % Flag indicating iteration limit has been reached
iter = 0;       % Number of iterations
x = x0;         % Initialize query point
while (threshFlag && iterFlag)
    delF = calcGradient(meritFunc, x, dx);
    
    % Calculate parameter increment
    diffF = delF - delF0;
    for itr = 1:num_params
%         gamma(itr) = dx(itr)*diffF(itr);
        gamma(itr) = gamma(itr)*diffF(itr); % DEBUG
    end
    gamma = abs(gamma)/norm(diffF)^2;
    
    x = x - gamma.*delF; % Update query point
    delF0 = delF;        % Update previous gradient
    
    % Check stopping conditions
    iter = iter+1;
    if (iter > iterLim) % Check iteration limit
        iterFlag = 0;
        exitFlag = 1;
    end
    % Check threshold limit
    if ( (norm(delF)/abs(F0) < thresh) && (iter > 1) )
        threshFlag = 0;
        exitFlag = 0;
    end
    
end