function [lambdak,V,flag] = JacobiRot(Ak,e,imax)
%JacobiRot execute the Jacobi rotation

flag = true;
if Transp(Ak) - Ak ~= 0
    flag = false;
    lambdak = [];
    V = [];
    return
end
error = inf;
count = 0;
len = length(Ak);
V = eye(len);
p = 0;
q = 0;
U = eye(len);

while and(error > e,count < imax)
    if count > 0
        Ak = Ak1; 
    end
    Amax(Ak);
    phi = (Ak(q,q) - Ak(p,p))/(2*Ak(p,q));
    
    if phi
        if phi < 0
            signal = -1;
        else
            signal = 1;
        end
        t = 1/(phi+signal*sqrt(phi^2+1));
    else
        t = 1;
    end
    c = 1/sqrt(1+t^2);
    s = t/sqrt(1+t^2);
    Amax(Ak);
    
    JacobiU;
    V = MatrixMulti(V,U);
    Ak1 = Transp(U)*Ak*U;
    
    Amax(Ak1);
    error = abs(Ak1(p,q));
    count = count + 1;
end
lambdak = zeros(1,len);
for i=1:len
    lambdak(i) = Ak1(i,i);
end

% Nested functions
    function Amax(A)
    amax = 0;
    for i=1:len
        for j=1:len
            if and(i ~= j,abs(A(i,j)) > amax)
                amax = A(i,j);
                p = i;
                q = j;
            end
        end
    end
    end

    function JacobiU
    U = eye(len);
    U(p,p) = c;
    U(q,q) = c;
    U(p,q) = s;
    U(q,p) = -s;
    end
end

