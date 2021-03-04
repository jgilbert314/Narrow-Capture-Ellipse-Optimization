function [ minQuan ] = meritFuncGeneral(x, InputStruct)

% Miscellaneous parameters
epsilon = InputStruct.epsilon;
D = InputStruct.D;
N = InputStruct.N;

% Domain specific parameters
a = InputStruct.a; b = InputStruct.b;
f = InputStruct.f;
beta = InputStruct.beta; xi_b = InputStruct.xi_b;
omega = InputStruct.omega;

% Value assigned when input violates problem constraints
% Do not use Inf. Program will eventually crash
badVal = 1e3;

% For some reason the optimizer generates NaN inputs. Don't know why.
if (any( isnan(x) ))
    x( isnan(x) ) = 0;
end


% Convert input elliptical coords to cartesian coords
if (a == b)
    [xC, yC] = calcPolar2Cartesian(x(1:N), x(N+1:end));
else
    [xC, yC] = calcElliptical2Cartesian(x(1:N), InputStruct.mu, x(N+1:end));
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

% Check if traps are in contact (shouldn't be necessary when checking for
% AMFPT <= 0)
if (epsilon > 0)
    for itrA = 1:N-1
        xA = x_vec(itrA, :);
        for itrB = itrA+1:N
            xB = x_vec(itrB, :);
            dist = sqrt(sum( (xA - xB).^2 ));
            if (dist <= 2*epsilon)
                disp('Traps in contact');
                disp(num2str(dist));
                minQuan = badVal;
                return;
            end
        end
    end
end

% Interaction matrix and coefficients
if (a == b) % Circular case
    G_mat = greensMatFuncCircle(x_vec, N);
else
    G_mat = greensMatFunc(x_vec, a, b, f, omega, beta, xi_b);
end

if (epsilon > 0)
    A = coeffVecFunc(G_mat, N, omega, epsilon, D);
else
    A = coeffVecFunc([], N, omega, epsilon, D);
end

% Quantity to be minimized
minQuan = sum( G_mat*A );


% % Physical solutions are those which have AMFPT > 0, as this value is a
% % time
% const = -1/D/N/log(epsilon); % Constant term used to calculate AMFPT
% AMFPT = const + minQuan;
% if (AMFPT <= 0)
%     disp('Bad AMFPT') % DEBUG
%     minQuan = badVal;
%     return;
% end


end