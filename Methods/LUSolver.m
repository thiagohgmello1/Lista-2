function [x,det,flag] = LUSolver(A,b)
%LUSolver Solve a linear system using LU decomposition

flag = true;

[L,U,P,det,flag2] = LUDecomp(A);
if flag2 == false
    flag = false;
    return
end

Pb = MatrixMulti(P,b);
y = LSDiagInf(L,Pb);
x = LSDiagSup(U,y);

end

