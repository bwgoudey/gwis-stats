function [logBM,   logPhi1, logPhi2] = logBM4curve( table1, table2,p1s, p2s, x1, x2)

logPhi1 = logP(table1,'logLow', p1s, x1)  ;
logPhi2 =   logP(table2,'logUp', p2s, x2);
logBM = max(logPhi1 + logPhi2);
end


function logPhi = logP(table, LowUp, p, x)

[ v, ix] = ismember(x, table.xVals);

%%% find flanking pVals
p2 = [ table.pVals,    table.pVals*0+1;...
    p,              p*0];
[V, ord] = sort(p2(:, 1));
p2 = p2(ord,:);
ip1 = cumsum(p2(:, 2));
I = find(p2(:,2)==0);

ip1 = max(1, ip1(I));   %  lower boundary
ip2 = min(ip1+1, length(table.pVals) );

t = table.(LowUp);

%%% allocate an average value:
logPhi = t(ip1, ix) .*( table.pVals(ip2)  -p)  + ...
        t(ip2, ix) .*( p - table.pVals(ip1)  );
logPhi = logPhi ./ max(0.000001, table.pVals(ip2) - table.pVals(ip1));
lopgPhi(t(ip1, ix) == -Inf |  t(ip2, ix) == -Inf) = -Inf;
end

