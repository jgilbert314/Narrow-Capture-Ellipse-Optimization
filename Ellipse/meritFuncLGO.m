function [ minQuan ] = meritFuncLGO(x_in)
% Author: Jason Gilbert
% Date: June 17, 2021
% Last Updated: 
% Version: 
%
% Summary:
%   Debugging function for use with LGO

% User Input
eps_0 = 0.025;
ecc = 0.125; % Eccentricity


N = length(x_in)/2; % Number of traps

% Vector of trap sizes (two traps twice as large as the rest)
epsilon = ones(N, 1);
epsilon(1:N-2) = eps_0;
epsilon(N-1:end) = eps_0*2;
nu = log(epsilon);

D = 1; % Coefficient of diffusivity

% Domain specific parameters
[ InputStruct ] = buildInputStruct(ecc, eps_0);
a = InputStruct.a; b = InputStruct.b;
f = InputStruct.f;
beta = InputStruct.beta; xi_b = InputStruct.xi_b;
omega = InputStruct.omega;

% Check if traps are common size
commSize_flag = (size(nu) == 1);
commSize_flag = commSize_flag(1);

% Value assigned when input violates problem constraints
% Do not use Inf. Program will eventually crash
badVal = 1e3;

% For some reason the optimizer generates NaN inputs. Don't know why.
if (any( isnan(x_in) ))
    x_in( isnan(x_in) ) = 0;
end


% Convert input elliptical coords to cartesian coords
if (a == b)
    [xC, yC] = calcPolar2Cartesian(x_in(1:N), x_in(N+1:end));
else
    [xC, yC] = calcElliptical2Cartesian(x_in(1:N), InputStruct.mu, x_in(N+1:end));
end

% Convert to row vector
[xR, ~] = size(xC);
if (xR ~= 1)
    xC = xC';
end
[yR, ~] = size(yC);
if (yR ~= 1)
    yC = yC';
end
% Convert input from vector to array
x_vec = [ xC' , yC' ];

% % Check if trap is within domain (shouldn't be necessary when using
% % elliptical coords)
% bound_check = (x_vec(:, 1)/a).^2 + (x_vec(:, 2)/b).^2; % Definition of ellipse boundary (must be <1)
% if ( any(bound_check > 1) )
%     minQuan = badVal;
%     return;
% end

% Interaction matrix and coefficients
if (a == b) % Circular case
    G_mat = greensMatFuncCircle(x_vec, N);
else        % Elliptical case
    G_mat = greensMatFunc(x_vec, a, b, f, omega, beta, xi_b);
end

if ( (commSize_flag) && all(nu == 0) )
    A = coeffVecFunc([], N, omega, nu, D);
else
    A = coeffVecFunc(G_mat, N, omega, nu, D);
end

% Quantity to be minimized
minQuan = sum( G_mat*A );


% AMFPT_approx = omega/4/pi^2/D/min(nu);
% if ( (minQuan + AMFPT_approx) <  AMFPT_approx*0.4 )
%     disp('Non-physical AMFPT') % DEBUG
%     minQuan = badVal;
%     return;
% end


% Check if traps are in contact (shouldn't be necessary when checking for
% AMFPT <= 0)
dist_tol = 2; % Factor of trap radii defining seperation
if (epsilon > 0)
    for itrA = 1:N-1
        xA = x_vec(itrA, :);
        radCoord = (xA(1)/a)^2 + (xA(2)/b)^2;
        
        % Check if trap contacts boundary
%         if (boundDist <= epsilon(itrA)*dist_tol^2)
        if (radCoord >= 0.9) % TESTING
            disp('Trap contacting boundary');
            minQuan = badVal;
            return;
        end
        
        % Check if traps are in contact with each other
        for itrB = itrA+1:N
            xB = x_vec(itrB, :);
            dist = sqrt(sum( (xA - xB).^2 ));
            if ( dist <= (epsilon(itrA) + epsilon(itrB))*dist_tol )
                disp('Traps in contact');
                disp(num2str(dist));
                minQuan = badVal;
                return;
            end
        end
    end
end


% Physical solutions are those which have AMFPT > 0, as this value is a
% time
% AMFPT = calcAMFPT(omega, D, N, G_mat, nu, A);
if ( minQuan <= -1 )
    disp('Non-physical AMFPT') % DEBUG
    minQuan = badVal;
    %%% TESTING
    fileID = fopen('TestLog.log', 'a');
    fprintf(fileID, '%0.15f,', [xC, yC]);
    fprintf(fileID, '\n');
    fclose(fileID);
    %%% END TESTING
    return;
end


end