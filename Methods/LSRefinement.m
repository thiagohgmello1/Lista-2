function [x,count,e,flag] = LSRefinement(A,x0,b,error, maxiter)
%LSRefinement refine a LS with x0 been the first approximation

r = b - MatrixMulti(A,x0);
c = LUSolver(A,r);
x = x0 + c;
e = Linfty(x,x0);
count = 0;

while (e > error) && (count < maxiter)
    x0 = x;
    r = b - MatrixMulti(A,x0);
    c = LUSolver(A,r);
    x = x0 + c;
    e = Linfty(x,x0);
    count = count + 1;
end
if count == maxiter
    flag = false;
    return;
end

end

