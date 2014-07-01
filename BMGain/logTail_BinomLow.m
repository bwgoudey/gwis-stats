function logtail   = logTail_BinomLow(x,t,p)
% computes lower tail for the binomial distribution
% tail = sum_{i=0}^x{n \choose i} p^i(1-p)^{t-i}
% used procedure in Appendix of "Filtering Genome Wide Interactions: 
% Classification Approach", D:\AdamPapers\2012\GWIS_Filtering/GWIS_filter.tex
%written: A.Kowalczyk, 14 June, 2012


EPS = 0.00001;
if p<= EPS
    logtail = 0;
    return
elseif p> 1-EPS
    logtail = -Inf ;
    return
end

%%% find arg x of maximal term
xdag =   ceil(t*p);
if  xdag*(1-p)/p/(t-xdag+1) >1
    xdag = xdag-1;
end
xdag = min(x, xdag);

persistent xdag0 x0 fMi0 fPl0
if isempty(xdag0) || xdag~=xdag0 || x~=x0
%if 1
    i = [1 : xdag];
    fMi0 =  (xdag-i+1)./(t-xdag+i);
    
    ii = [  1 : x - xdag];
    fPl0 =  (t-xdag-ii+1)./(xdag+ii);
    
    x0 = x;
    xdag0 = xdag;
end
    A = p/(1-p);

if xdag == 0
    fMi = 0;
else    
    fMi = fMi0 /A ;
    %fMi = [1; fMi];  % append entry for f_0
    fMi = cumprod(fMi);
end
if xdag == x
    fPl = 0;
else
    fPl = fPl0 *A ;
    fPl = cumprod(fPl);
end
 logtail = logBinom(t,xdag,p) + log(1 + sum(fMi) + sum(fPl));
end

