function [I3] = Simpson38(h,data)
%Simpson38 calculate integral using 3/8 Simpson's rule integration

len = length(data);
I3 = 0;

for j=1:len
    i = j-1;
    
    if (i == 0 || j == len)
        c = 1;
    elseif (~rem(i,3))
        c = 2;
    else
        c = 3;
    end
    I3 = I3 + c*data(j);
end

I3 = I3*3*h/8;

end

