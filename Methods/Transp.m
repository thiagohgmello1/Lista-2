function [At, flag] = Transp(A)
%Transp determine transposition at A matrix

flag = true;
A_size = size(A);

% Check if transposition is possible (A is a square matrix)
if A_size(1) ~= A_size(2)
    flag = false;
    At = [];
    return;
end

len = size(A);
At = zeros(len);

for i=1:len
    At(i,:) = A(:,i);
    At(:,i) = A(i,:);
end

end

