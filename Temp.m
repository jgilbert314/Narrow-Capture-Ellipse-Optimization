clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input Variables
a = 1; b = 0.5;         % Dimensions of ellipse
x = a/2; y = b/2;       % Evaluation point coordinates
x_0 = a/3; y_0 = b/3;   % Reference point coordinates

% Vector of trap coordinates
x_vec = [ ...
    0,  a/2 ; ...
    0, -a/3 ; ...
     b/3, 0 ; ...
    -b/3, 0 ...
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Domain dependent terms
omega = pi*a*b;
% Equation (4.3b)
f = sqrt(a^2 - b^2);
beta = (a - b)/(a + b);
xi_b = -log(beta)/2;

% Micellaneous terms
[N_0, ~] = size(x_vec);
epsilon = [1/20, 1/20, 1/20, 1/20]';
D = 1;
% coeff = omega/2/pi/D/N; % DEBUG - uncomment
coeff = omega/2/pi/D; % DEBUG - delete
nu = -1./log(epsilon);

I_mat = eye(N_0);             % Identity matrix

% Calculate AMFPT
G_mat = greensMatFunc(x_vec, a, b, f, omega, beta, xi_b);



%%

% Define right-hand side of equation
RHS = nu;

% Define left-hand side of equation
LHS = zeros(N_0);
for itrA = 1:N_0
    for itrB = 1:N_0
        
        term = 2*pi*nu(itrA)*G_mat(itrA, itrB);
        if (itrA == itrB)
            term = term + 1;
        end
        LHS(itrA, itrB) = term;
    
    end
end

% Calculate vector of coefficients
A_prop = LHS\RHS;
corr_factor = sum(A_prop)/coeff;
A = A_prop/corr_factor;

A = coeffVecFunc(G_mat, N_0, omega, nu, D);

%%

% Calculate AMFPT
AMFPT = coeff/nu + 2*pi/N_0 * ones(1, N_0)*G_mat*A;

minQun = sum( G_mat*A );
A_unit = A/norm(A);


%% Ref
nu_ref = nu(1);
A_ref = coeffVecFunc(G_mat, N_0, omega, nu_ref, D);
AMFPT_ref = coeff/nu_ref + 2*pi/N_0 * ones(1, N_0)*G_mat*A;

minQun_ref = sum( G_mat*A_ref );
A_unit_ref = A_ref/norm(A_ref);

plot(A - A_ref)

%% Calc AMFPT

% nu_mat = diag(nu);
% LHS = 2*pi*nu_mat*(G_mat*A) + A;
% nu./LHS

AMFPT = calcAMFPT(omega, D, N_0, G_mat, nu, A);
AMFPT_ref = calcAMFPT(omega, D, N_0, G_mat, nu_ref, A_ref);


%% Calc Trap Scalings

eps_0 = 1/20;
N_0 = 10;

N_p = N_0;
k = 1;

[ eps_p ] = calcTwoTrapSizes(eps_0, N_0, N_p, k);

eps_1 = sqrt( N_0/(N_0 + N_p*(k^2 - 1)) )*eps_0;
eps_2 = k*eps_1;

A_0 = N_0*pi*eps_0^2;
A_1 = (N_0 - N_p)*pi*eps_1^2;
A_2 = N_p*pi*eps_2^2;

A_0 - (A_1 + A_2)

A_p = (N_0 - N_p)*pi*eps_p(1)^2 + N_p*pi*eps_p(2)^2;
A_0 - A_p