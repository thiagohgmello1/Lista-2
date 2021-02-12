function [x,count,error,flag] = JacobiSolver(A,b,conv,e,kmax)
%JacobiSolver solve a LS using Jacobi iterative method

% Square matrix evaluation
A_size = size(A);
if A_size(1) ~= A_size(2)
    x = [];
    flag = false;
    return
end
len = A_size(1);

% Variables initialization
xk = zeros(len,1);
J = zeros(A_size);
c = zeros(len,1);
count = 0;
error = inf;

% J, x0 and c matrices determination
for i=1:len
    for j=1:len
        if i ~= j
            J(i,j) = -A(i,j)/A(i,i);
        end
    end
    xk(i) = b(i)/A(i,i);
    c(i) = b(i)/A(i,i);
end

% Convergence garantee method
if or(conv == "n", conv == "N")
    ratio = SpecRatio(J,"PowerMet",e,kmax);
    if abs(ratio) > 1
        flag = false;
        x = [];
        return
    end
elseif or(conv == "s", conv == "S")
    flag = DiagDom(A);
    if flag == false
        x = [];
        return
    end
else
    flag = false;
    return
end

% Jacobi method
while and(error >= e, count < kmax)
    xk1 = MatrixMulti(J,xk) + c;
    
    count = count + 1;
    error = Linfty(xk1,xk)/abs(max(xk));
    xk = xk1;
end
x = xk1;

end

