function [I1] = TrapComp(h,data)
%TrapComp calculate integral using trapezoidal integration

len = length(data);

c = ones(1,len);
c = c*2;
c(1) = 1;
c(len) = 1;

I1 = sum(c.*data)*h/2;

end

