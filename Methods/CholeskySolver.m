function [x, flag] = CholeskySolver(A,b)
%CholeskySolver Solve a linear system using Cholesky decomposition

flag = true;

[L,Lt,flag2] = CholeskyDecomp(A);
if flag2 == false
    flag = false;
    return
end

y = LSDiagInf(L,b);
x = LSDiagSup(Lt,y);

end

