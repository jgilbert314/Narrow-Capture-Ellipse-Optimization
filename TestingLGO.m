clear;

addpath('C:\LGO');
addpath('Ellipse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Input
N = 10;
epsilon = 0.025; % Trap sizes (hardcoded in meritFuncLGO)
ecc = 0.125; % Eccentricity of domain (hardcoded in meritFuncLGO)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LGO Options

% Maximum allowed value (limited by 32-bit program?)
% 1e9 too big
max_val = 1e4;
% Search type
% (0, 1, 2, 3) -> (local, branch+bound, global adaptive random, multi-start)
optimizationType = 3;
% Max number of function evaluations
% Not sure how any of these work. Not well documented
globalSample = max_val; localSample = max_val; numberTrying = max_val;

constraintWeight = 1; % Not used, required argument
min_fun_val = 0; % Threshold for solution to be considered acceptable?
tol_val = 0;     % Local search improvement target (not sure how it works)
randomSeed = 0;  % Seed for random search
searchTime = max_val; % Maximum allowable search time

options = [optimizationType, globalSample, localSample, numberTrying ...
    , constraintWeight, min_fun_val, min_fun_val, tol_val, 1e-6, 1e-6, randomSeed, searchTime];

%% Run LGO

nvars = 2*N;
% Generate initial coordinates
% Polar coordinates of spiral
theta = linspace(0, N, N);
theta = mod(theta, 2*pi); % Wrap polar coordinates
[ InputStruct ] = buildInputStruct(ecc, epsilon);
r = linspace(0, 0.999*InputStruct.A, N);
x0 = [theta, r];
% Actual bounds enforced within merit function
lbounds = zeros(1, nvars);
ubounds = [2*pi*ones(1, N) , ones(1, N)];
bnds = [lbounds ; x0 ; ubounds]';
ncons = 0; ctypes = [];

[fval, opt_param, ~, addInfo] = matlablgo('meritFuncLGO', nvars, bnds, ncons, ctypes, options);