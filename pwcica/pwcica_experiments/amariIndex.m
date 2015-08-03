function ai = amariIndex(A,B)

scalingA = repmat(sqrt(mean(inv(A).^2))',[1 size(A,2)]);
scalingB = repmat(sqrt(mean(inv(B).^2))',[1 size(B,2)]);

a = (scalingA.*A)/(scalingB.*B);
M = size(A,1);
ai = 1/(2*M*(M-1))*sum(sum(abs(a),2)./max(abs(a'))'-1) ...
    + 1/(2*M*(M-1))*sum(sum(abs(a))./max(abs(a))-1);

end