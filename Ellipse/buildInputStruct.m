function [ InputStruct ] = buildInputStruct(ecc, epsilon)

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

[A, mu] = calcEllipticalParameters(a, b);

% Input Structure setup
InputStruct = struct();

% Miscellaneous parameters
InputStruct.epsilon = epsilon;
InputStruct.D = D;
InputStruct.nu = nu;

% Domain specific parameters
InputStruct.a = a; InputStruct.b = b;
InputStruct.f = f; 
InputStruct.beta = beta; InputStruct.xi_b = xi_b;
InputStruct.omega = omega;

% Elliptical Coords
InputStruct.A = A;
InputStruct.mu = mu;

end