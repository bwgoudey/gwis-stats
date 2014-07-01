function ct=getCTfromPlink(prb, bed_file, label)
    %Wrapper around extractPlinkGenos and geno2Ct 

    nSample=length(label);
    nGeno=3;
    nPheno=length(unique(label));
    assert(nPheno==2, 'unable to handle nPheno != 2');
    
    max_prb_per_iter=10000;
    show_pbar=size(prb, 1)>max_prb_per_iter;
    if(show_pbar)
        h=textprogressbar(0, 'getCTfromPlink: ');    
    end
   
    ct=nan(size(prb,1), nPheno * nGeno^size(prb,2));
    for prb_start = 1:max_prb_per_iter:size(prb, 1)        
        if(show_pbar)
            h=textprogressbar(prb_start/size(prb, 1), h);
        end        
        curr_I = prb_start:min(prb_start+max_prb_per_iter-1, size(prb, 1));
        
        genos=extractPlinkGenos(bed_file, prb(curr_I,:), nSample, 'unpacked');
        ct(curr_I, :)=genos2Ct(prb(curr_I,:), genos, label, nGeno);        
    end
    if(length(prb)>max_prb_per_iter)
        h=textprogressbar(1, h);
        h=textprogressbar(sprintf('DONE\n'), h);
    end    
    assert(sum(isnan(ct(:)))==0);
end