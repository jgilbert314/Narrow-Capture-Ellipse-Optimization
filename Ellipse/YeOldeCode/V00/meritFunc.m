function [ minQuan ] = meritFunc(x, InputStruct)

% Convert x to row vector
[xR, ~] = size(x);
if (xR ~= 1)
    x = x';
end


% tic;

% Miscellaneous parameters
epsilon = InputStruct.epsilon;
D = InputStruct.D;
N = InputStruct.N;

% Domain specific parameters
a = InputStruct.a; b = InputStruct.b;
f = InputStruct.f; 
beta = InputStruct.beta; xi_b = InputStruct.xi_b;
omega = InputStruct.omega;

% Convert input from vector to array
x_vec = [ x(1:N)' , x(N+1:end)' ];

% Check if trap is within domain
bound_check = (x_vec(:, 1)/a).^2 + (x_vec(:, 2)/b).^2; % Definition of ellipse boundary (must be <1)
if ( any(bound_check > 1) )
    minQuan = Inf;
    return;
end

% Interaction matrix and coefficients
G_mat = greensMatFunc(x_vec, a, b, f, omega, beta, xi_b);
A = coeffVecFunc(G_mat, N, omega, epsilon, D);

% Quantity to be minimized
minQuan = sum( G_mat*A );

% CompTime = toc;

% % Print user feedback
% disp(['Inputs: ', num2str(x)]);
% disp(['Merit:  ', num2str(minQuan)]);
% disp(['Computation Time: ', num2str(CompTime)]);
% 
% fprintf(InputStruct.fileID, '%0.15f, ', x);
% fprintf(InputStruct.fileID, '%0.15f, ', minQuan);
% fprintf(InputStruct.fileID, '%0.15f\n', CompTime);


end