function [ fltGSS, sel,fltPP ] = calcMxGSS( ct )

    minCover=0.02;
            

    [~, nV] = size(ct);
    nV = nV/2;
    assert( nV==round(nV));

%     prevPermut = zeros(size(ct,1),9);
% 
%     % preallocate
%     sel.fltGSS_cntr = zeros(size(ct,1), 1);
%     sel.fltGSS_prtv = zeros(size(ct,1), 1);
%     sel.OR_cntr = zeros(size(ct,1), 1);
%     sel.sen_cntr = zeros(size(ct,1), 1);
%     sel.x0_cntr = zeros(size(ct,1), 1);
%     sel.x1_cntr = zeros(size(ct,1), 1);            
%     sel.OR_prtv = zeros(size(ct,1), 1);
%     sel.spe_prtv = zeros(size(ct,1), 1);
%     sel.x0_prtv = zeros(size(ct,1), 1);
%     sel.x1_prtv = zeros(size(ct,1), 1);
%     sel.num_cntr = zeros(size(ct,1), 1);
%     sel.num_prtv = zeros(size(ct,1), 1);
%     sel.prevPermut = zeros(size(ct,1),1);

%     if(all(ct(ii,:) == 0))
%         ct(ii,1)    = sum(ct(1,1:nV));
%         ct(ii,nV+1) = sum(ct(1,1+nV:end));
%     end



    [ct0, ct1]=seperateCtrlCasesCt(ct);
    t0 = sum( ct0,2 ) ;   % this takes into account missing values
    t1 =  sum( ct1,2 ) ;
    
    
    mxGSS('primeGSS', max(t0), max(t1), 1000, minCover);
    [cntrIdx, gsstab] = mxGSS('gss', ct0', ct1');
    cntrIdx = cntrIdx';
    gsstab = -gsstab';
%%

    I_cntr = repmat(1:9, size(gsstab, 1), 1);
    tmpIdx = repmat(cntrIdx, 1, 9);
    I_cntr = I_cntr >= tmpIdx;
    [sel.fltGSS_cntr, i_cntr] = max(gsstab .* I_cntr, [], 2);
    [sel.fltGSS_prtv, i_prtv] = max(gsstab .* ~I_cntr, [], 2);

    sel.prevPermut=getPrevPerm(ct0, ct1);
    prevIdx = prevPerm2Idx(sel.prevPermut);
    x0 = cumsum(applyPerm(ct0, prevIdx, 2), 2);
    x1 = cumsum(applyPerm(ct1, prevIdx, 2), 2);

    x0_cntr = x0(i_cntr);
    x1_cntr = x1(i_cntr);

    sel.num_cntr = i_cntr;
    sel.num_prtv = 9 - i_prtv;

    sel.OR_cntr = x1_cntr./(t1-x1_cntr).*(t0 - x0_cntr)./x0_cntr;
    sel.sen_cntr = x1_cntr./t1;
    sel.x0_cntr = x0_cntr;
    sel.x1_cntr = x1_cntr;

    x0_prtv = x0(i_prtv);
    x1_prtv = x1(i_prtv);
    sel.x0_prtv = x0_prtv;
    sel.x1_prtv = x1_prtv;

    sel.OR_prtv =  1./( x1_prtv./(t1-x1_prtv).*(t0 - x0_prtv)./x0_prtv);
    sel.spe_prtv = 1 - x0_prtv./t0;
    fltPP = gsstab;
   
    %histSS = [];
    fltGSS = max([sel.fltGSS_cntr, sel.fltGSS_prtv],[], 2);
    sel.fltGSS = fltGSS;
end

