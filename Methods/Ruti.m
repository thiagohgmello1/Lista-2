function [lambdak,flag] = Ruti(Ak,e,imax)
%Ruti execute the Rutishauser method

error = inf;
count = 0;
len = length(Ak);

while and(error > e,count < imax)
    [L,U] = LUDecomp(Ak);
    Ak = MatrixMulti(U,L);
    MaxError;
    count = count + 1;
end

lambdak = zeros(1,len);
for i=1:len
    lambdak(i) = Ak(i,i);
end

if count > 50
    flag = false;
    return
else
    flag = true;
end

    function MaxError
        maxerror = 0;
        for i=1:len
            for j=1:i
                if and(abs(Ak(i,j)) > maxerror, i ~= j)
                    maxerror = Ak(i,j);
                end
            end
        end
        error = maxerror;
    end

end