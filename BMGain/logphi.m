function val = logphi(p1,p2,x1,x2,t1,t2)
%
% computes the natural log of phi
% Written: A.Kowalczyk, 15 june, 2012

val   = logTail_BinomLow(x1,t1,p1) + logTail_BinomUp(x2,t2,p2);

