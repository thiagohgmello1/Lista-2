function [error, flag] = Linfty(v,w)
%Finfty determine the infinity norm btween two arrays

flag = true;
lenv = length(v);
lenw = length(w);

if lenv ~= lenw
    error = 0;
    flag = false;
    return
end
error = 0;

for i=1:lenv
    e = abs(v(i) - w(i));
    if e > error
        error = e;
    end
end

end

