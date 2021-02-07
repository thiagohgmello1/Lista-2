function [q,Ainv,detA,flag] = PCLevFad(A)
%PCLevFad Determine characteristic polynomial, matrix inverse and
%determinant

flag = true;
A_size = size(A);
if A_size(1) ~= A_size(2)
    Ainv = [];
    q = [];
    detA = 0;
    flag = false;
    return;
end

len = A_size(1);
Bn = eye(len);
I = eye(len);
q = zeros(1,(len + 1));
q(1) = 1;

for i=1:len
    qn = 0;
    Bn_1 = Bn;
    An = MatrixMulti(A,Bn_1);
    for k=1:len
        qn = qn + An(k,k);
    end
    q(i + 1) = -qn/i;
    Bn = An - (qn/i)*I;
end

if Bn ~= zeros(len)
    flag = false;
    Ainv = [];
    detA = 0;
    return;
end

Ainv = (1/q(len + 1))*Bn_1;
detA = q(len + 1)*(-1)^len;

end

