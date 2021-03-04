function [ x0 ] = generateInitialConfiguration( N )
    
    % Polar coordinates of spiral
    theta = linspace(0, N, N);
    r = linspace(0.1, 0.75, N);
    
    x = r.*cos(theta);
    y = r.*sin(theta);
    
    x0 = [x, y];
    
end