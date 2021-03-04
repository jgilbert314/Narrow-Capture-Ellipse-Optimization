function [ r, theta ] = calcCartesian2Polar(x, y)

r = sqrt(x.^2 + y.^2);
theta = atan2(y, x);

end