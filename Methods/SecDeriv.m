function [deriv2] = SecDeriv(data)
%SecDeriv defines second order derivative at point xk using tree points.
%
% The entrances are data(1) == f(x1), data(2) == f(x2), data(3) == 
% f(x3), data(4) == h. x3 > x2 > x1.

y1 = data(1);
y2 = data(2);
y3 = data(3);
h = data(4);

deriv2 = (y3-2*y2+y1)/h^2;

end

