% clear;

nu = InputStruct.nu;
epsilon = InputStruct.epsilon;
D = InputStruct.D;

% Domain specific parameters
a = InputStruct.a; b = InputStruct.b;
f = InputStruct.f;
beta = InputStruct.beta; xi_b = InputStruct.xi_b;
omega = InputStruct.omega;

data_arr = readmatrix('TestLog.log');
data_arr = data_arr(:, 1:end-1);
[nR, nC] = size(data_arr);
N = nC/2;
colour_list = hsv(N);
colour_list_rows = hsv(nR);

% Boundary coords
bX = linspace(-a, a, 1e3);
bY = b*sqrt( 1 - (bX/a).^2 );


%% Convert coords
coord_arr = zeros(nR, 2*N);
for itr = 1:nR
    
    x = data_arr(itr, 1:N);
    y = data_arr(itr, N+1:end);
    if (a == b)
        [xC, yC] = calcPolar2Cartesian(x, y);
    else
        [xC, yC] = calcElliptical2Cartesian(x, InputStruct.mu, y);
    end
    
    coord_arr(itr, 1:N) = xC;
    coord_arr(itr, N+1:end) = yC;
    
end

%% Mess
bound_arr = zeros(nR, N);

% nR = 10;
minR = zeros(1, nR);
for itr = 1:nR
    x = data_arr(itr, 1:N);
    y = data_arr(itr, N+1:end);
    % Convert input elliptical coords to cartesian coords
    if (a == b)
        [xC, yC] = calcPolar2Cartesian(x, y);
    else
        [xC, yC] = calcElliptical2Cartesian(x, InputStruct.mu, y);
    end
    x_vec = [xC ; yC]';
    
    bound_arr(itr, :) = ( x_vec(:, 1).^2/a^2 + x_vec(:, 2).^2/b^2 );
    
    hold on
%     for itrN = 1:N
%         plot(xC(itrN), yC(itrN), '.', 'MarkerEdgeColor', colour_list(itrN, :));
%     end
%     plot(bound_arr(itr, :), '.');
%     drawnow;
    
    thisMin = inf;
    for itrA = 1:N
        for itrB = itrA+1:N
            thisR = sqrt(sum( (x_vec(itrA, :) - x_vec(itrB, :)).^2 ));
            if (thisR < thisMin)
                thisMin = thisR;
            end
        end
    end
    minR(itr) = thisMin;
%     bound = ( x_vec(:, 1).^2/a^2 + x_vec(:, 2).^2/b^2 );
%     if (any(bound >= 0.90))
%         disp('Close to bound');
%     end
   

    
end

%% Plot coords (bad)
hold off
for itrN = 1:N
  
    if (itr == 2)
        hold on
    end
    plot(data_arr(:, itrN), data_arr(:, N+itrN), '.', 'MarkerEdgeColor', colour_list(itrN, :));
    
    drawnow;
end
hold off

%% Plot Coords
hold off
plot(bX, bY, 'k');
hold on
plot(bX, -bY, 'k');
axis equal
for itr = 1:nR
    plot(data_arr(itr, 1:N), data_arr(itr, N+1:end), '.', 'MarkerEdgeColor', colour_list_rows(itr, :));  
end

%% Calc coeff vec

A_arr = zeros(nR, N);
for itr = 1:nR
    x = data_arr(itr, 1:N);
    y = data_arr(itr, N+1:end);
    x_vec = [x', y'];
    G_mat = greensMatFunc(x_vec, a, b, f, omega, beta, xi_b);
    A_arr(itr, :) = coeffVecFunc(G_mat, N, omega, nu, D);
end



%% Bound check

xC = data_arr(:, 1:N);
yC = data_arr(:, N+1:end);

bound_check = (xC/a).^2 + (yC/b).^2;
hold on
for itr = 1:nR
    plot( bound_check(itr, :), '.' );
end

%%

% Boundary coords
a_p = 0.9*a; b_p = 0.9*b;
bX_p = linspace(-a_p, a_p, 1e3);
bY_p = b_p*sqrt( 1 - (bX_p/a_p).^2 );

hold on
plot(bX_p, bY_p, 'r');
plot(bX_p, -bY_p, 'r');

%% Check merit output

InputStruct.N = N;
% Trap size parameters
epsilon_0 = InputStruct.epsilon_0;
N_p = InputStruct.N_p;
k = InputStruct.k;
InputStruct.epsilon = calcTwoTrapSizes(epsilon_0, N, N_p, k); % TESTING

merit_vals = zeros(1, nR);
for itr = 1:nR
    merit_vals(itr) = meritFuncGeneral(data_arr(itr, :), InputStruct);
end

%%

data = [0.222617977172110,0.767361813192162,0.783342142976210,0.551601026852761,-0.770929085682921,-0.251455637268565,0.906182679183581,0.681501512929176,0.684030770027528,0.393306927435040,-0.057900562307050,0.273373832701862,-0.259075725538851,-0.709215497156406,0.102239020883126,0.726091572250533,-0.142697335717978,-0.490088186579391,-0.024310983129436,-0.804761819838286];

x = data(1:N);
y = data(N+1:end);

plot(x, y, '.');
hold on
plot(bX, bY, 'k');
plot(bX, -bY, 'k');
hold off

%%

InputStruct.N = N;
% Trap size parameters
epsilon_0 = InputStruct.epsilon_0;
N_p = InputStruct.N_p;
k = InputStruct.k;
InputStruct.epsilon = calcTwoTrapSizes(epsilon_0, N, N_p, k); % TESTING
InputStruct.nu = -1./log(InputStruct.epsilon);

meritFuncGeneral(data, InputStruct);

%%

% Boundary coords
bX = linspace(-a, a, 1e3);
bY = b*sqrt( 1 - (bX/a).^2 );

hold on
plot(bX, bY, 'k');
plot(bX, -bY, 'k');
plot(xC, yC, 'k.')

test = (xC/a).^2 + (yC/b).^2;
max(test)

%%

hold on
for itr = 1:nR
   x = data_arr(itr, 1:N);
   y = data_arr(itr, N+1:end);
   rad = (x/a).^2 + (y/b).^2;
   rad = sort(rad);
   
   plot(rad);
end
hold off