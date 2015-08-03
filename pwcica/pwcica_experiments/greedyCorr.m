function [out,revInd] = greedyCorr(data1,data2)

m = size(data1,1);
n = size(data2,1);
% N = size(data1,2);

cM = abs(corr(data1',data2'));
out = zeros(1,m);
revInd = zeros(1,m);
for ii = 1:m
    [best,indx] = max(cM(:));
    [rr,cc] = ind2sub(size(cM),indx);
    cM(rr,:) = zeros(1,n);
    cM(:,cc) = zeros(m,1);
    out(rr) = best;
    revInd(rr) = cc;
end
    




end