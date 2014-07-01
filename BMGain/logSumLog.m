function logSum = logSumLog(log1, log2)
% computes sum of two numbers with given logs
%{
    logSum  = NaN * log1;
    
    II = log1 == log2;
    logSum(II) = log1(II) + log(2);

    I = log1 > log2;
    logSum(I) = log1(I) + log1p(exp(log2(I)-log1(I)));

    II = log1 < log2;
    logSum(II) = log2(II) + log1p(exp(log1(II)-log2(II)));
%}
    if isvector(log1)
        vecShape = size(log1);
        C = max([log1(:), log2(:)], [], 2); % always reshape to use columns         
        logSum = C + log1p(exp(-abs(log1 - log2)));
        II = ~isfinite(C);
        logSum(II) = C(II) + log(2);
        logSum = reshape(logSum, vecShape);
    else
        logSum  = NaN * log1;
    
        II = log1 == log2;
        logSum(II) = log1(II) + log(2);

        I = log1 > log2;
        logSum(I) = log1(I) + log1p(exp(log2(I)-log1(I)));

        II = log1 < log2;
        logSum(II) = log2(II) + log1p(exp(log1(II)-log2(II)));
    end
end
 
