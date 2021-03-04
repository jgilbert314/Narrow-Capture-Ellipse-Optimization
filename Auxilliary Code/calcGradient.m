function [ delF ] = calcGradient(meritFunc, x0, dx)

num_params = length(x0);
x = x0;

F0 = meritFunc(x0);
F = zeros(1, num_params);

for itr = 1:num_params
    x(itr) = x(itr) + dx(itr);
    F(itr) = meritFunc(x);
    x(itr) = x0(itr);
end

delF = (F - F0)./dx;

end