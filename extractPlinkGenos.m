function geno=extractPlinkGenos(bed_file, prb, nSample, return_format)
% Call:
%    geno=getPlinkGenos(file, snp_indicies, nSample)
% Description:
%    Read in genotypes of provided list of probes and convert into 
%    additive encoding (ie {0,1,2} indicating occurrence of minor allele)
%  Input:
%    bed_file = path to bed file
%    prb = indicies of probes to get genotype data for
%    nSample = number of samples in the given dataset
%  Output:
%    geno structure containing
%      .prb = list of probe indicies in descending order
%      .data = genotype data (size=[prb, nSample]) corresponding to each prb
    if(~exist('return_format', 'var'))
        return_format='packed';
    end
    %% 
    assert(~isempty(prb), 'List of probes is empty');
    prb=unique(reshape(prb, numel(prb),1));
    
    %err='provided vector should be of individual probes only';
    %assert(size(prb,2)==1, err);
    
    err=sprintf('provided bed file does not exist:\n\t%s\n', bed_file);
    assert(exist(bed_file, 'file')~=0, err);
    
    geno_per_byte=4;
    magic_bytes=3; %Part of plink specification
    m=memmapfile(bed_file, 'Offset', magic_bytes);
    
    %I figured it may be easier to read in data if probes are sorted
    prb=sort(prb);

    %Convert index range into byte range
    bytes_per_snp=ceil(nSample/geno_per_byte);
    geno_data=[];
    geno.prb=[];
    max_prb_per_iter=10000;
    
    if(length(prb)>max_prb_per_iter)
        h=textprogressbar(0, 'Extract Genos: ');    
    end
   
    for prb_start = 1:max_prb_per_iter:length(prb)        
        if(length(prb)>max_prb_per_iter)
            h=textprogressbar(prb_start/length(prb), h);
        end
        curr_prb = prb(prb_start:min(prb_start+max_prb_per_iter-1, length(prb)));
        byte_indicies_start=((curr_prb-1)*bytes_per_snp)+1;
        byte_indicies_end=byte_indicies_start+bytes_per_snp-1;
        %Return range between each start and end.
        byte_indicies=cell2mat(arrayfun(@(x,y) (x:y),...
            byte_indicies_start',byte_indicies_end', 'UniformOutput', 0));

        %Read in all data from bed file
        geno_data=[geno_data; m.Data(byte_indicies')];
        geno.prb=[geno.prb; curr_prb];
    end
    
    %Make sure we got our iteration correct
    assert( sum(geno.prb~=prb)==0, 'iteration through loop was incorrect')    
    
    if(length(prb)>max_prb_per_iter)
        h=textprogressbar(1, h);
        h=textprogressbar(sprintf('DONE\n\r'), h);
    end
    
    if(strcmp(return_format, 'packed'))
        geno.data=uint8(geno_data);
    elseif(strcmp(return_format, 'unpacked'))
        geno.data=convertToAdditiveEncoding(geno_data,nSample,bytes_per_snp,geno_per_byte);
    end
end


function data=convertToAdditiveEncoding(geno_data,nSample,bytes_per_snp,geno_per_byte)
% Convert to encoding based on the number of times the minor allele occurs 
% [0,1,2]. Missing value in unpacked data 3. 
    plink_mapping=[2,3,1,0];
    %Create list of all possible plink 8 byte values
    geno_combs=combn(plink_mapping,length(plink_mapping));
    geno_combs=geno_combs(:, length(plink_mapping):-1:1);
    %uint16 and +1 so geno_data can act as a index in 1-based array
    geno_data=uint16(geno_data);
    geno_data=reshape(geno_combs(geno_data+1, :)',...
        bytes_per_snp*geno_per_byte, [])';

    %Remove extra 0s at end from padding to whole byte
    data=geno_data(:, 1:nSample);
end
