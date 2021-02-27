function [deriv] = deriv2(data)
% deriv2 determines the derivative at xk of a f function using two points.
%
% The entrances are data(1) == f(x1), data(2) == f(x2), data(3) == h.
% x2 > x1. 
% If dir == "P" or "p", the derivative is in the progressive direction.
% Else if dir == "R" or "r", the derivative is in the regressive direction

y1 = data(1);
y2 = data(2);
h = data(3);

deriv = (y1 - y2)/h;

end

