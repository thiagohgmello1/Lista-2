function [lambdak1,yk1] = PowerMet(A,e,imax)
%PowerMet determine the higher eigenvalue of A matrix

count = 0;
A_size = size(A);
len = A_size(1);
lambdak = zeros(len,1);
lambdak1 = lambdak;
yk = randi(10,1,len);
zk1 = MatrixMulti(A,yk);
a = max(abs(zk1));
yk1 = 1/a*zk1;
error = inf*ones(len,1);

while and(error >= e, count < imax)
    zk1 = MatrixMulti(A,yk);
    a = max(abs(zk1));
    yk1 = 1/a*zk1;
end


end

