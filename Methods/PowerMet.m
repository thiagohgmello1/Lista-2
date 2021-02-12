function [lambda,yk] = PowerMet(A,e,imax)
%PowerMet determine the higher eigenvalue of A matrix

count = 0;
A_size = size(A);
len = A_size(1);
error = inf;

yk = Transp(randi(10,1,len));
while all(yk) == false
    yk = Transp(randi(10,1,len));
end

zk = MatrixMulti(A,yk);
a = max(abs(zk));
yk = 1/a*zk;
zk = MatrixMulti(A,yk);
lambdak1 = zk./yk;

while and(error >= e, count < imax)
    lambdak = lambdak1;
    a = max(abs(zk));
    yk = 1/a*zk;
    zk = MatrixMulti(A,yk);
    lambdak1 = zk./yk;
    error = abs(max((lambdak-lambdak1)./lambdak1));
    yk = (1/a)*zk;
    count = count + 1;
end
lambda = mean(lambdak1);

end

