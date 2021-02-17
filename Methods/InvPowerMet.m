function [lambda,yk] = InvPowerMet(A,e,imax)
%PowerMet determine the higher eigenvalue of A matrix

count = 0;
A_size = size(A);
len = A_size(1);
error = inf;

yk = Transp(randi(10,1,len));
while all(yk) == false
    yk = Transp(randi(10,1,len));
end

[L,U,P] = LUDecomp(A);

Py = MatrixMulti(P,yk);
y = LSDiagInf(L,Py);
zk = LSDiagSup(U,y);

a = max(abs(zk));
yk = 1/a*zk;

Py = MatrixMulti(P,yk);
y = LSDiagInf(L,Py);
zk = LSDiagSup(U,y);
a = max(abs(zk));
yk = 1/a*zk;

lambdak1 = zk./yk;

while and(error >= e, count < imax)
    lambdak = lambdak1;
    
    Py = MatrixMulti(P,yk);
    y = LSDiagInf(L,Py);
    zk = LSDiagSup(U,y);
    
    lambdak1 = zk./yk;
    error = abs(max((lambdak-lambdak1)./lambdak1));
    yk = (1/a)*zk;
    count = count + 1;
end
lambda = 1/mean(lambdak1);

end

