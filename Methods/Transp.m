function [At] = Transp(A)
%Transp determine transposition at A matrix

A_size = size(A);

At = zeros(A_size(2),A_size(1));

for i=1:A_size(1)
    At(:,i) = A(i,:);
end

for j=1:A_size(2)
    At(j,:) = A(:,j);
end

end

