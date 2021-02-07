function [C, flag] = MatrixMulti(A,B)
% MatrixMulti realize multiplication between two matrix.

% Initializaitons
flag = true;
Asize = size(A);
Bsize = size(B);
C = zeros(Asize(1),Bsize(2));
aux = 0;

% Check if operation is possible
if Asize(2) ~= Bsize(1)
    flag = false;
    return;
end

% Matrix multiplication
for i=1:Asize(1)
    for j=1:Bsize(2)
        for k=1:Asize(2)
            aux = aux + A(i,k)*B(k,j);
        end
        C(i,j) = aux;
        aux = 0;
    end
end
end

