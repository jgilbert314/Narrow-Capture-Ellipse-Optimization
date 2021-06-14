clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input Variables
a = 1; b = 0.5;         % Dimensions of ellipse
x = a/2; y = b/2;       % Evaluation point coordinates
x_0 = a/3; y_0 = b/3;   % Reference point coordinates

% Vector of trap coordinates
x_vec = [ ...
    0,  a/2 ; ...
    0, -a/3 ...
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Domain dependent terms
omega = pi*a*b;
% Equation (4.3b)
f = sqrt(a^2 - b^2);
beta = (a - b)/(a + b);
xi_b = -log(beta)/2;

% Micellaneous terms
[N, ~] = size(x_vec);
epsilon = 0.05;
D = 1;
coeff = omega/2/pi/D/N;
nu = -1/log(epsilon);

% Calculate AMFPT
G_mat = greensMatFunc(x_vec, a, b, f, omega, beta, xi_b);
A = coeffVecFunc(G_mat, N, omega, epsilon, D);
AMFPT = coeff/nu + 2*pi/N * ones(1, N)*G_mat*A;

minQun = sum( G_mat*A );
A_unit = A/norm(A);

%% Initialization
clear;

% Input variables
ecc = 0.3;                % Eccentricity of ellipse
N = 25;                   % Number of traps
epsilon = 0;              % Trap size
filename = 'logP.matopt'; % Logfile name

% Domain dependent terms
[a, b] = calcAxis(ecc); % Axis of ellipse
omega = pi*a*b;
% Equation (4.3b)
f = sqrt(a^2 - b^2);
beta = (a - b)/(a + b);
xi_b = -log(beta)/2;

% Micellaneous terms
D = 1;
coeff = omega/2/pi/D/N;
nu = -1/log(epsilon);


%% Minimize AMFPT
InputStruct = struct();

% Miscellaneous parameters
InputStruct.epsilon = epsilon;
InputStruct.D = D;
InputStruct.N = N;

% Domain specific parameters
InputStruct.a = a; InputStruct.b = b;
InputStruct.f = f; 
InputStruct.beta = beta; InputStruct.xi_b = xi_b;
InputStruct.omega = omega;

% PSWARM Options
opts = optiset('maxiter', 1e3, 'maxfeval', 1e9, ...
    'tolrfun', 1e-3, 'tolafun', Inf, ...
    'display', 'iter');
opts.solverOpts = pswarmset('swarm_size', 100);

% % fminsearch Options
% opts = optimset('fminsearch');
% opts.MaxFunEvals = 1e4;
% opts.MaxIter = 1e4;

% x0 = [ linspace(-a*0.5, a*0.5, N), linspace(-b*0.5, b*0.5, N) ]'; % Initial positions
x0 = [ (rand(1, N) - 0.5)*a, (rand(1, N) - 0.5)*b]';
B = [a*ones(1, N), b*ones(1, N)]';
LB = -B; UB = B;   % Parameter bounds (wrong, need to be linear constraints)

OptHandle = tic;
InputStruct.fileID = fopen(filename, 'w+');
try
    mFunc = @(x)meritFunc(x, InputStruct);
%     [opt_param, fval, exitflag, output] = fmincon(mFunc, x0, [], [], [], [], LB, UB);
%     [opt_param, fval, exitflag, output] = fminsearch(mFunc, x0, opts);
    [opt_param, fval, exitflag, output] = opti_pswarm(mFunc, LB, UB, x0, [], [], opts);
catch ME
    fclose(InputStruct.fileID);
    rethrow(ME);
end
fclose(InputStruct.fileID);
CompTime = toc(OptHandle);


% Plotting
res = 100;
x_plt = linspace(-a, a, res);
plt_quan = sqrt( b^2*(1 - x_plt.^2/a^2) );

x0_plt = x0;
xT_plt = opt_param;

hold off
plot(x_plt, plt_quan, 'k');
hold on
plot(x_plt, -plt_quan, 'k');
plot(xT_plt(1:N), xT_plt(N+1:end), 'b .'); % Final position
plot(x0_plt(1:N), x0_plt(N+1:end), 'k *'); % Initial position
hold off

% ax_lim = 1.1;
% xlim([-a*ax_lim, a*ax_lim]);
% ylim([-b*ax_lim, b*ax_lim]);
axis equal;
