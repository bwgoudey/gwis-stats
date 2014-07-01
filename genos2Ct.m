function ct=genos2Ct(prb, geno, y, nGeno)
% genos2Ct
% CALL
%     ct=genos2Ct(prb, geno, y, nGeno)
% DESCRIPTION
%     takes genotypes in an additive-encoding (prob. from extractPlinkGenos)
%     and covert into contingency tables for the specified probes
%
% SEE
%   extractPlinkGenos
%
    assert(nnz(prb)/size(prb,2)==size(prb,1), 'not all probes are >0');
    %Covert missing to majority
    geno.data(geno.data==3)=NaN;
    
    nGenoCombs=power(nGeno,size(prb, 2));
    ct=zeros(size(prb, 1), nGenoCombs*length(unique(y)));
    
    max_prb_per_iter=1e9;
    show_pbar=size(prb, 1)>max_prb_per_iter;
    if(show_pbar)
        h=textprogressbar(0, 'Count Genos: ');    
    end
   
    for prb_start = 1:max_prb_per_iter:size(prb, 1)        
        if(show_pbar)
            h=textprogressbar(prb_start/size(prb, 1), h);
        end        
        curr_I = prb_start:min(prb_start+max_prb_per_iter-1, size(prb, 1));

        %Convert SNps to their place in the genotype data
        [a, geno_ind]=ismember(prb(curr_I, :), geno.prb);
        assert(sum(sum(geno_ind<1))==0)
    
        %if 2 genos are passed in, then 
        switch size(prb, 2)
            case 1
                g=geno.data(geno_ind, :);
                %ct(:, :, 1)=histc(g(y==1), 0:(nGenoCombs-1));
                %ct(:, :, 2)=histc(g(y==2), 0:(nGenoCombs-1));
            case 2
                g=geno.data(geno_ind(:, 2), :)*nGeno+geno.data(geno_ind(:, 1), :);
            case 3
                g1=geno.data(geno_ind(:, 1), :);
                g2=geno.data(geno_ind(:, 2), :);
                g3=geno.data(geno_ind(:, 3), :);
                g=g3*nGeno^2+g2*nGeno+g1;
            case 4
                g1=geno.data(geno_ind(:, 1), :);
                g2=geno.data(geno_ind(:, 2), :);
                g3=geno.data(geno_ind(:, 3), :);
                g4=geno.data(geno_ind(:, 4), :);
                g=g4*nGeno^3+g3*nGeno^2+g2*nGeno+g1;                
        end
        %%
    
        ct(curr_I, 1:nGenoCombs)      =histc(g(:, y==1), 0:(nGenoCombs-1), 2);
        ct(curr_I, (nGenoCombs+1):end)=histc(g(:, y==2), 0:(nGenoCombs-1), 2);
        %assert(nnz(sum(ct(curr_I, 1:9),2)==0)==0, 'all CT entries are missing')
    end

    if(show_pbar)
        h=textprogressbar(1, h);
        h=textprogressbar(sprintf('DONE\n'), h);
    end        
end
