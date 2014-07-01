function  [logBM, x1s, x2s]   = logBM_table( table1, table2,p1s, p2s, x1s, x2s )

if isempty(x1s)
    x1s = table1.xVals;
end
if isempty(x2s)
    x2s = table2.xVals;
end

ip1 = findIndexLow(p1s, table1 );
%ip1 = min(max(table1.pVals), ip1+1);  % thi sis more conservative

ip2 = findIndexLow(p2s, table2 );
ip2 = min(length(table2.pVals), ip2+1);  % thi sis more conservative

[TF1, ix1 ] = ismember(x1s, table1.xVals);
[TF2, ix2 ] = ismember(x2s, table2.xVals);
if ~all(TF1 &TF1)
    error('logBM_table: missing values of x')
end

table1.logLow = table1.logLow(ip1,ix1);
table2.logUp = table2.logUp(ip2,ix2);

logBM =  zeros(length(x1s), length(x2s));

for i1 = 1 : length(ix1)
    M = repmat(table1.logLow(:, i1), 1, length(ix2));
    M = M + table2.logUp;
    logBM(i1, :) = max(M,[],1);
end

end

function ip = findIndexLow(p, table)

%%% find flanking pVals
p2 = [ table.pVals,    table.pVals*0+1;...
    p,              p*0];
[V, ord] = sort(p2(:, 1));
p2 = p2(ord,:);
ip1 = cumsum(p2(:, 2));
I = find(p2(:,2)==0);

ip = max(1, ip1(I));   %  lower boundary
end