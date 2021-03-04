function [ G ] = calcGreensFunc(x, y, x_0, y_0, a, b)

r = sqrt(x^2 + y^2);
r_0 = sqrt(x_0^2 + y_0^2);

% Equation (4.3b)
f = sqrt(a^2 - b^2);
beta = (a - b)/(a + b);
xi_b = -log(beta)/2;

xi = xiFunc(x, y, f);         % Equation (4.4a)
xi_0 = xiFunc(x_0, y_0, f);
eta = etaFunc(x, y, f);       % Equation (4.4b)
eta_0 = etaFunc(x_0, y_0, f);

z_arr = z_arrFunc(xi, xi_0, xi_b, eta, eta_0); % Equation (4.4b)

% Greens function
xi_ang = max( [xi, xi_0] );
G = greensFunc(r, r_0, a, b, beta, xi_ang, z_arr);

end