function [L,Lt,flag] = CholeskyDecomp(A)
%CholeskDecomp define Cholesky decomposition from A matrix (A = LL')

[At,flag] = Transp(A);

if or(ne(At,A),(flag ~= true))
    L=[];
    flag = false;
    return
end

len = length(A);

% First column calculation
L(1,1) = sqrt(A(1,1));
for i=2:len
    L(i,1) = A(i,1)/L(1,1);
end

Lsum = 0;
for j=2:len
    for k=1:(j-1)
        Lsum = Lsum + L(j,k)^2;
    end
    L(j,j) = sqrt(A(j,j) - Lsum);
    Lsum = 0;
    for i=(j + 1):len
        for k=1:(j-1)
            Lsum = Lsum + L(i,k)*L(j,k);
        end
        L(i,j) = (A(i,j) - Lsum) / L(j,j);
        Lsum = 0;
    end
end
Lt = Transp(L);

end

