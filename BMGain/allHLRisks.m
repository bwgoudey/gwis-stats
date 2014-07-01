function [P0, P1] = allHLRisks( p0s, p1s)


k = length(p0s);
%%%generate binary expansion
B = [0;1];
for ic = 2 : k
    B = [B, zeros(size(B(:,1))); B, ones(size(B(:,1)))];
end
P0a = B;
P1a = B;
for ic=1 : k
   P0a(:,  ic) = P0a(:,  ic)  *  p0s(ic);
   P1a(:,  ic) = P1a(:,  ic)  *  p1s(ic);
end
P0 = sum(P0a,2);
P1 = sum(P1a,2);
end