

function sel = runGssAnalysis(plink_files, snp_list, GSS_type)

%% Setup variables
if(~exist('GSS_type','var'))
    GSS_type='GSS';
end
%% Load in and convert data
snp_info=readBim([plink_files '.bim'], '\t');
snp_info.prb=(1:length(snp_info.snp_id))';
samp_info=readFam([plink_files '.fam'], ' ');

%% Filter if need be
if(exist('snp_list', 'var'))
    if(size(snp_list,2)==1)
        nUniqSnp = length(unique(snp_list));
        total_pairs = nchoosek(nUniqSnp,2);
    elseif(size(snp_list,2)==2)
        nUniqSnp = length(unique(snp_list));
        total_pairs = size(snp_list,1);
    else
        assert(0, 'list poorly specified');
    end
else
    nUniqSnp=size(snp_info.prb,1);
    total_pairs = nchoosek(nUniqSnp,2);
end

%% Create list of selected pairs (sel) 
%  Assume that all pairs are currently being recorded. 
if(strcmp(GSS_type, 'MxGSS'))
    step_size = 9000;
else
    step_size = 10;
end

prev_pair = [1,1];

sel.prb=zeros(total_pairs,2);
sel.GSS=zeros(total_pairs,1);
sel.snp_id=cell(total_pairs,2);

%%
h=textprogressbar(0, sprintf('Calc GSS (%d): ', total_pairs));

for i = 1:step_size:total_pairs   
    %%
    h=textprogressbar(i/total_pairs, h);    
    curr_step = (i:min(i+step_size-1, total_pairs))';
    
    %Figure out the prb index and the snp identifier that we are looking at
    if(size(snp_list,2)==2)
        [~, prb]=ismember(snp_list(curr_step, :), snp_info.snp_id);   
        sel.snp_id(curr_step,:)=snp_list(curr_step, :);
    elseif(size(snp_list,2)==1)
        pairs = convertToPairInd(prev_pair,step_size, nUniqSnp);
        snp_pairs = snp_list(pairs);
        [~, prb]=ismember(snp_pairs, snp_info.snp_id);        
        sel.snp_id(curr_step,:)=snp_pairs;
    else    
        pairs = convertToPairInd(prev_pair,step_size, nUniqSnp);
        prb = snp_info.prb(pairs);
        sel.snp_id(curr_step,:)=snp_info.snp_id(pairs);
    end
    sel.prb(curr_step,:) = prb;
    
    % Get contingency table and calculated GSS
    ct9x2 = getCTfromPlink(prb, [plink_files '.bed'], samp_info.label);   
        
    if(strcmp(GSS_type, 'MxGSS'))
        sel.GSS(curr_step,:) = calcMxGSS(ct9x2);
    else
        sel.GSS(curr_step,:) = calcGSS(ct9x2);
    end
    
    if(size(snp_list,2)~=2)
        prev_pair = pairs(end,:);
    end
end
textprogressbar(sprintf('DONE\n'),h);
%%
[~,base]=fileparts(plink_files);

save([base '_GSS.mat'], 'sel');
% 
% fid=fopen([base '_GSS.txt'], 'w');
% for i = 1:length(sel.GSS)
%     rs_1=snp_info.snp_id{sel.prb(i,1)};
%     rs_2=snp_info.snp_id{sel.prb(i,2)};
%     fprintf(fid, '%s\t%s\t%d\t%d\t%.4f\n', rs_1,rs_2,sel.prb(i,:), sel.GSS(i));    
% end
% fclose(fid);
end
