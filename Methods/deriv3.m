function [deriv] = deriv3(dir,data)
%deriv3 defines first order derivative at point xk using tree points.
%
% The entrances are data(1) == f(x1), data(2) == f(x2), data(3) == 
% f(x3), data(4) == h. x3 > x2 > x1.

deriv = inf;

y1 = data(1);
y2 = data(2);
y3 = data(3);
h = data(4);

if (dir == "R" || dir == "r")
    deriv = 1/(2*h)*(-3*y1+4*y2-y3);
elseif (dir == "C" || dir == "c")
    deriv = 1/(2*h)*(y3-y1);
elseif (dir == "P" || dir == "p")
    deriv = 1/(2*h)*(y1-4*y2+3*y3);
end

end

