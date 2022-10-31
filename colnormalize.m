function [ B ] = colnormalize ( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[m,n] = size(A);
B = zeros(m,n);
if n>1
    for k=1:n
        meanA = mean(A(:,k));
        stdA = std(A(:,k));
        B(:,k) = (A(:,k)-meanA)/stdA;
    end
else
    meanA = mean(A);
    stdA = std(A);
    B = (A-meanA)/stdA;
end
end

