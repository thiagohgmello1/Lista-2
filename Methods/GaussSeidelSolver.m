function [xk1,count,error,flag] = GaussSeidelSolver(A,b,e,kmax)
%GaussSeidelSolver solve a LS using Gauss-Seidel iterative method

flag = true;
% Square matrix evaluation
A_size = size(A);
if A_size(1) ~= A_size(2)
    xk1 = [];
    flag = false;
    return
end
len = A_size(1);

% Variables initialization
D = zeros(len);
E = zeros(len);
F = zeros(len);

xk = zeros(len,1);
xk1 = zeros(len,1);
count = 0;
error = inf;

% First xk determination
for i=1:len
    xk(i) = b(i)/A(i,i);
end

% D, E and F matrices
for i=1:len
    for j=1:len
        if i == j
            D(i,i) = A(i,i);
        elseif i < j
            E(i,j) = -A(i,j);
        elseif i > j
            F(i,j) = -A(i,j);
        end
    end
end

% S and d matrices
[q,DEinv] = PCLevFad(D-E);
S = MatrixMulti(DEinv,F);
d = MatrixMulti(DEinv,b);

% Gauss-Seidel method
while and(error >= e, count < kmax)
    xk1 = MatrixMulti(S,xk) + d;
    count = count + 1;
    error = Linfty(xk1,xk)/abs(max(xk));
    xk = xk1;
end

end

