function [ fltGSS, sel,sel_all,fltPP ] = calcGSS( ct )

    minCover=0.02;


    [~, nV] = size(ct);
    nV = nV/2;
    assert( nV==round(nV));

    prevPermut = zeros(size(ct,1),9);

    % preallocate
    sel.fltGSS_cntr = zeros(size(ct,1), 1);
    sel.fltGSS_prtv = zeros(size(ct,1), 1);
    sel.OR_cntr = zeros(size(ct,1), 1);
    sel.sen_cntr = zeros(size(ct,1), 1);
    sel.x0_cntr = zeros(size(ct,1), 1);
    sel.x1_cntr = zeros(size(ct,1), 1);            
    sel.OR_prtv = zeros(size(ct,1), 1);
    sel.spe_prtv = zeros(size(ct,1), 1);
    sel.x0_prtv = zeros(size(ct,1), 1);
    sel.x1_prtv = zeros(size(ct,1), 1);
    sel.num_cntr = zeros(size(ct,1), 1);
    sel.num_prtv = zeros(size(ct,1), 1);
    sel.prevPermut = zeros(size(ct,1),1);

    nCalcs=size(ct,1);
    max_iter_for_pbar=  1e5;  
    if(nCalcs>max_iter_for_pbar)
        h=textprogressbar(0, sprintf('Calc GSS (%d): ',size(ct,1)));
    end
    for ii = 1 : size(ct,1);
        if(nCalcs>max_iter_for_pbar)
            h=textprogressbar(ii/size(ct,1), h);
        end
        if(all(ct(ii,:) == 0))
            ct(ii,1)    = sum(ct(1,1:nV));
            ct(ii,nV+1) = sum(ct(1,1+nV:end));
        end

        t0 = sum( ct(ii,1:nV ),2 ) ;   % this takes into account missing values
        t1 =  sum( ct(ii, nV + [1:nV]),2 ) ;

        ct0 = reshape( ct(ii, 1:9), 3, 3 );
        ct1 = reshape( ct(ii, 10:18), 3, 3 );

        [lgBM, snp_ab, x0, x1, rRisk] = gainBM( ct0, ct1, [], 0.005, [num2str(ii) ': '],'binomDirect');    

        sel.prevPermut(ii,:) =getPrevPerm(ct(ii, 1:9), ct(ii, 10:18));

        prevPermut(ii,:) = snp_ab(:,3)';
        sen = x1/t1;
        spe = 1-x0/t0;

        %%% eliminate too small cover
        lgBM(spe < minCover | sen< minCover) = 0;

        I_cntr = spe >= sen;
        I_prtv = sen >= spe;

        %bwgoudey 1/10/13: Pretty certain this (-lgBM + 1)
        %inflates the GSS p-value by an order of magnitude. 
        [sel.fltGSS_cntr( ii,1 ), i_cntr]  = max((-lgBM) .* I_cntr );  % "+1" is to avoid issue with 0 value
        [sel.fltGSS_prtv( ii,1 ), i_prtv]  = max((-lgBM) .* I_prtv );

        x0_cntr = x0(i_cntr);
        x1_cntr = x1(i_cntr);

        sel.num_cntr(ii,1) = i_cntr;
        sel.num_prtv(ii,1) = 10 - i_prtv;  % AK changed on 9Feb

        sel.OR_cntr( ii,1 ) = x1_cntr/(t1-x1_cntr)*(t0 - x0_cntr) /x0_cntr;
        sel.sen_cntr( ii,1 ) = x1_cntr/t1;
        sel.x0_cntr( ii,1 ) =x0_cntr;
        sel.x1_cntr( ii,1 ) = x1_cntr;

        x0_prtv = x0(i_prtv);
        x1_prtv = x1(i_prtv);
        sel.x0_prtv( ii,1 ) =x0_prtv;
        sel.x1_prtv( ii,1 ) = x1_prtv;

        sel.OR_prtv( ii,1 ) =  1/( x1_prtv/(t1-x1_prtv)*(t0 - x0_prtv) /x0_prtv);  % AK changed on 9 Feb, 2013
        sel.spe_prtv( ii,1 ) = 1 - x0_prtv/t0;

        fltPP(ii,:) = -lgBM';   % this is a vector of all values, for distribution creation


        sel_all=struct('fltGSS', lgBM, 'x0', x0, 'x1', x1, 'g_a', snp_ab(:, 1), 'g_b', snp_ab(:, 2), 'prevPermut', snp_ab(:,3));
    end
    if(nCalcs>max_iter_for_pbar)
        textprogressbar(sprintf('DONE\n'),h);
    end
    %histSS = [];
    fltGSS = max([sel.fltGSS_cntr, sel.fltGSS_prtv],[], 2);
    sel.fltGSS = fltGSS;
end

