function [I2] = Simpson13(h,data)
%Simpson1 calculate integral using 1/3 Simpson's rule integration

len = length(data);
I2 = 0;

for j=1:len
    i = j-1;
    
    if (i == 0 || j == len)
        c = 1;
    elseif (rem(i,2))
        c = 4;
    else
        c = 2;
    end
    I2 = I2 + c*data(j);
end

I2 = I2*h/3;

end

