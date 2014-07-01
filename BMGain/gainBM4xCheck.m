function [log10BM, snp_ab, x0, x1] = gainBM4xCheck( ct0, ct1, actions, dp, titlePrfx,figH )

%AK, 17 Sep, 2012, ouptput changed to log10BM and titlePrfx added

% actions = {'plot', 'printTable' 'plotCloud'}
% dp = 0.005; default
ct0 = double(ct0);
ct1 = double(ct1);

t0 = sum(sum(ct0));     %   control
t1 = sum(sum(ct1));     %   case
log_10 = log(10);

persistent t0_0 t1_0 dp0 ps table0 table1 SNP
if isempty(dp0) |   dp0 ~=dp  | t0 ~= t0_0 | t1 ~= t1_0
  %  ps = [0:0.001:0.1 0.1:0.005: 0.2 0.2:0.01:0.8 0.8:0.005: 0.9 0.9:0.001:1.0]';
    ps = [0:0.0001:0.002, 0.002 :0.001:0.1, 0.1:0.005: 0.2, 0.2:0.01:0.8, ...
          0.8:0.005: 0.9, 0.9:0.001:0.98, 0.98:0.0001:1.0]';
    ps = unique(ps);
    table0  = tailsTable('binomDirect' ,'Low',ps, [], t0, [],1);
    table1  = tailsTable('binomDirect' ,'Up' ,ps, [], t1, [],1);
    dp0 = dp; t0_0 = t0; t1_0 = t1;

    SNP= [  0     1     2     0     1     2     0     1     2;...
            0     0     0     1     1     1     2     2     2]';
end

%%% find snp_a
cta(:,1) = sum(ct0,2);
cta(:,2) = sum(ct1,2);
cta(:,3) = [0:2]';      % store genotype value

[va, ord]  = sort(cta(:,2)./ max(0.00001,cta(:,1))* t0/t1, 'descend');  
cta = cta(ord,:);
cta(:,4) = cumsum(cta(:,1));
cta(:,5) = cumsum(cta(:,2));

% col=1 - cnt0
% col=2 - cnt1
% col=3 - genotype
% col=4 - cumulated cnt0
% col=5 - cumulated cnt1
 
%%% find snp_b
ctb(:,1) = sum(ct0,1)';
ctb(:,2) = sum(ct1,1)';
ctb(:,3) = [0:2]';      % store genotype value

[vb, ord]  = sort(ctb(:,2)./ max(0.00001,ctb(:,1))* t0/t1, 'descend'); 
ctb = ctb(ord,:);
ctb(:,4) = cumsum(ctb(:,1));
ctb(:,5) = cumsum(ctb(:,2));

%%% find convex hull
% col=6 - snp 1==a, 2 ==b
% col=7 - diffCnt0
% col=8 - diffCnt1
% col=9 - difftangent

ct = [  cta, ones(3,1),    cta(:,4)/t0, cta(:,5)/t1 ; ...
        ctb, 2*ones(3,1),  ctb(:,4)/t0, ctb(:,5)/t1];
ct(:,9) = ct(:,8)./ max(0.0001, ct(:,7))  ;
ct(7,:) =  [ t0, t1, -1, t0, t1, -1, 1, 1, 1] ;
hall = [0,0];
while ~isempty(ct)
     [v, row] = max(ct(:,9));
     
     %%% add to the hall
     hall = [hall; ct(row,[ 4 5] )]; 
     %%% find reminder
     ct(:, 7) = (ct(:, 4)  - ct(row, 4))/t0;
     ct(:, 8) = (ct(:, 5)  - ct(row, 5))/t1;
     ct(:,9) = ct(:,8)./ max(0.0001, ct(:,7));
     ct = ct(ct(:,7) > 0, :);

end

%%% show both snps 
ct_ab = [reshape(ct0, 9,1), reshape(ct1, 9,1), [1:9]'];
risk_ab  =  ct_ab(:,2)./ ct_ab(:,1)* t0/t1;
[ v, ord] = sort(risk_ab, 'descend');

ct_ab = ct_ab(ord,:);
ct_ab0 = ct_ab;
ct_ab(:,1) = cumsum(ct_ab(:,1));
ct_ab(:,2) = cumsum(ct_ab(:,2));

%%% find boundary
param.H0Type = 'pieceLin';
param.p1 = hall(:,1)/t0;
param.p2 =  hall(:,2)/t1 ;
%param.title = ['H0Type= ' H0Type  ];
NW_p0 = [0: dp:1]';
NW_p1 =  NWboundary(NW_p0, param);

%%% find BM
for ii = 1 : 9 
    xx0 = ct_ab(ii,1);
    xx1 = ct_ab(ii,2);
    [logBM(ii,1), logPhi1, logPhi2] = logBM4curve( table0, table1, NW_p0,  NW_p1, xx0, xx1);    
end
 log10BM = logBM/log_10;

 x0=ct_ab(:, 1); 
 x1=ct_ab(:, 2);

 % patch for some obvious issues
I0 = find(x0==0 & x1==0 | x0==t0 & x1==t1);
logBM(I0,1) = 0;

[logBMst, iSt] = min(logBM);

 x0st=ct_ab(iSt,1); 
 x1st=ct_ab(iSt,2);

 snp_ab=[SNP( ct_ab(:,3),:), ct_ab(:,3)]; 
%% display
if ismember('plot', actions) 
    cta =   cta; % convert to %
    ctb =   ctb; % convert to %
    if nargin < 5  | isempty(figH)
     figH = figure(12); clf
    end
    hsp = subplot(1,1,1);
    plot(   [0;cta(:,4)/t0;1], [0;cta(:,5)/t1;1], 'b:v',...
            [0;ctb(:,4)/t0;1], [0;ctb(:,5)/t1;1], 'g:^',...
            [0;hall(:,1)/t0;1], [0;hall(:,2)/t1;1], 'm',...
            [0;ct_ab(:,1)/t0;1], [0;ct_ab(:,2)/t1;1], 'r-o',...
            'LineWidth',2,'MarkerSize', 5 ), hold on
    hl =    legend('snp_a', 'snp_b', 'convHall', 'snp_a & snp_b',...
                    'Location', 'SouthEast');
    if ismember('plotCloud', actions)
         [p0all_ab, p1all_ab] = allHLRisks(ct_ab0(:,1)/t0, ct_ab0(:,2)/t1);
         [p0all_a, p1all_a] = allHLRisks(cta(:,1)/t0, cta(:,2)/t1);
         [p0all_b, p1all_b] = allHLRisks(ctb(:,1)/t0, ctb(:,2)/t1);
         plot(p0all_ab, p1all_ab, 'ro', 'MarkerSize', 5)
         plot( p0all_a, p1all_a, 'bv', 'MarkerFaceColor', 'b','MarkerSize', 5);
         h=plot( p0all_b, p1all_b, 'g^', 'MarkerFaceColor', 'g','MarkerSize', 5);
     end
        %%% repeat lines for clarity
     plot(   [0;cta(:,4)/t0;1], [0;cta(:,5)/t1;1], 'b:v',...
            [0;ctb(:,4)/t0;1], [0;ctb(:,5)/t1;1], 'g:^',...
            [0;hall(:,1)/t0;1], [0;hall(:,2)/t1;1], 'm',...
            [0;ct_ab(:,1)/t0;1], [0;ct_ab(:,2)/t1;1], 'r-o',...
            'LineWidth',2,'MarkerSize', 5 ), hold on
     plot(  1-[0;ct_ab(:,1)/t0;1], 1-[0;ct_ab(:,2)/t1;1], 'r-',...
            1-[0;hall(:,1)/t0;1], 1-[0;hall(:,2)/t1;1], 'm',...
            [0,1], [0, 1], 'k--'), hold on
     axis square
%      xlabel(['controls (t_0=' num2str(t0) ')'])
%      ylabel(['cases (t_1=' num2str(t1) ')'])
     xlabel(['Low Risk (out of t_0=' num2str(t0) ' Controls)'])
     ylabel(['High Risk (out of t_1=' num2str(t1) ' Cases)'])
     
     OR = x1st/max(x0st, 0.00001)*(t0-x0st)/max(t1-x1st, 0.00001);
     MAF = x1st/t1;
     title([titlePrfx 'log_{10}BM=' num2str(logBMst/log(10), '%3.1f') ,...
            ', OR=' num2str(OR, '%3.1f') ', p_0^*=', num2str(100*x0st/t0, '%3.1f') '%, p_1^*=',...
            num2str(100*x1st/t1, '%3.1f') '%']);
     
     %%% mark iSt
     hp = plot(x0st/t0, x1st/t1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
     
     axis([0 1  0 1])
end
%% printTable
if ismember('printTable', actions)
    
    ORall = ct_ab(:,2)./max(ct_ab(:,1), 0.00001).*(t0-ct_ab(:,1))./max(t1-ct_ab(:,2)  , 0.00001);
    fprintf('\n  \tlog10(BM) \t OR \tx0 \tx1 \tp0%s \tp1%s \tsnp1 \tsnp2 \tsnp1&2', '%', '%');
    for ii = 1 : 9
        fprintf('\n%d \t%4.1f \t%4.2f \t%d \t%d \t%3.1f \t%3.1f \t%d \t%d \t%d', ...
            ii,  logBM(ii)/log_10, ORall(ii), x0(ii), x1(ii), ...
            100*x0(ii)/t0, 100*x1(ii)/t1,   snp_ab(ii,1),  snp_ab(ii,2), snp_ab(ii,3));
    end
    fprintf('\n')
end


end


