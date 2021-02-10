function [flag] = DiagDom(A)
%DiagDom determine if A is diagonal dominant

A_size = size(A);
if A_size(1) ~= A_size(2)
    flag = false;
    return
end
len = A_size(1);
Lsum = 0;

for i=1:len
    diag = abs(A(i,i));
    for j=1:len
        if j ~= i
            Lsum = Lsum + abs(A(i,j));
        end
    end
    if diag <= Lsum
        flag = false;
        return 
    end
    Lsum = 0;
end
flag = true;

end

