function [x,flag] = LSDiagInf(A,b)
%LSDiagSup Solve a superior diagonal linear system

A_size = size(A);
flag = true;

if A_size(1) ~= A_size(2)
    flag = false;
    return;
end
len = A_size(1);
xsum = 0;

x = zeros(len,1);

for i=1:len
    for j=1:i
        xsum = xsum - A(i,j)*x(j);
    end
    x(i) = (xsum + b(i))/A(i,i);
    xsum = 0;
end

end

