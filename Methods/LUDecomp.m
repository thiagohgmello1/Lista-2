function [L,U,P,det,flag] = LUDecomp(A)
% LUDecomp realize LU decomposition determining L, U and P matrix.

% Inicializations
flag = true;
A_size = size(A);

% Check if decomposition is possible (A is a square matrix)
if A_size(1) ~= A_size(2)
    flag = false;
    return;
end

% Initializations
det = 1;
detL = 1;
detU = 1;
count = 0;
L = zeros(A_size(1));
Laux = L;
U = A;
Uaux = U;
P = zeros(A_size(1));
len = A_size(1);
Psum = zeros(len,1);

for j=1:len
    
    pos = 0;
    pivot = 0;
    
    % P' matrix determination
    for i=1:len
        if and((abs(U(i,j)) > abs(pivot)),(Psum(i) == 0))
            pivot = U(i,j);
            pos = i;
        end
    end
    if pos == 0
        flag = false;
        break;
    end
    P(pos,j) = 1;
    Psum = sum(P,2);
    
    % L and U matrix determination
    for i=1:len
        if Psum(i) == 0
            m = -U(i,j) / pivot;
        else
            m = 0;
        end
        L(i,j) = -m;
        
        for j2=j:len
            if Psum(i) == 0
                U(i,j2) = U(i,j2) + m * U(pos,j2);
            end
        end
    end
    L(pos,j) = 1;
end

% Check if the system is determined
if flag == false
    return
end

% Build P, L and U matrix after pivonting process
for j=1:len
    for i=1:len
        if P(i,j) == 1
            auxU = U(i,:);
            auxL = L(i,:);
            Uaux(j,:) = auxU;
            Laux(j,:) = auxL;
        end
    end
end
U = Uaux;
L = Laux;
P = Transp(P);

% Determinant calculation
Paux = P;
for j=1:len
    for i=1:len
        if Paux(i,j) == 1 && i ~= j
            count = count + 1;
            auxP = Paux(j,:);
            Paux(j,:) = Paux(i,:);
            Paux(i,:) = auxP;
        end
    end
end

for i=1:len
    detL = detL * L(i,i);
    detU = detU * U(i,i);
end
det = detL*detU*(-1)^count;

end

