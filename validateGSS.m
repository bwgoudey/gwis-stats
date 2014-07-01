

cohorts={'UK1b', 'UK2b', 'FIN', 'NL', 'IT'};
proj_dir='/media/ntfs_2TB/Research/gwis-pp/projects/cel/';
getPlinkName = @(cohort) sprintf('%s/data/cel%s_noX_HH/A_raw/cel%s_noX_HH', proj_dir, cohort,cohort);

getBadGenoName = @(cohort) sprintf('%s/data/cel%s_noX_HH/B_annt/cel%s_noX_HH_bad_geno.mat', proj_dir, cohort,cohort);

%% Load data
plink_files = getPlinkName(cohorts{1});
snp_info=readBim([plink_files '.bim'], '\t');

%% Limit to chr 6    
chr=6;
min_pos=25e6;
max_pos=34e6;
tf_pos = snp_info.chr==chr & snp_info.pos >= min_pos & snp_info.pos <= max_pos;

%% Remove badly genotyped SNPs
rs_bad = {};
for i = 1:length(cohorts)   
    if(exist(getBadGenoName(cohorts{i}), 'file'))
        load(getBadGenoName(cohorts{i})); % ==> bad_geno;     
        rs = cellstr([repmat('rs',length(bad_geno.rs),1) num2str(bad_geno.rs)]);
        rs_bad = [rs_bad; regexprep(rs, '\s+', '')];        
    end
end
tf_good_geno=~ismember(snp_info.snp_id,rs_bad);

%% Run analysis 
sel_MxGSS = runGssAnalysis(plink_files, snp_info.snp_id(tf_pos&tf_good_geno), 'MxGSS');
signif_thresh = -calcBonf(length(snp_info.chr),2,0);
tf_signif = sel_MxGSS.GSS>signif_thresh-2;
sel = runGssAnalysis(plink_files, sel_MxGSS.snp_id(tf_signif,:), 'GSS');

sel=cutStruc2Inds(sel, find(sel.GSS>signif_thresh));

[~,base]=fileparts(plink_files);
save([base '_GSS.mat'], 'sel');

sels{1}=sel;
clear sel;

%% Get signif in other cohorts. 
signif = zeros(length(sels{1}.snp_id),length(cohorts));
signif(:,1) = ones(length(sels{1}.snp_id),1);

for i = 2:length(cohorts)
    plink_files = getPlinkName(cohorts{i});
    snp_info=readBim([plink_files '.bim'], '\t');
    
    tf = all(ismember(sels{1}.snp_id, snp_info.snp_id), 2);
    sels{i} = runGssAnalysis(plink_files, sels{1}.snp_id(tf,:), 'GSS');
    
    
    %Convert snp_ids to numbers
    all_snp_id = union(sels{1}.snp_id,sels{i}.snp_id) ;
    [~,ai] = ismember(sels{1}.snp_id,all_snp_id) ;
    [~,bi] = ismember(sels{i}.snp_id,all_snp_id) ;
    [tf,loc] = ismember(bi,ai,'rows');
    
    signif_thresh = -calcBonf(length(snp_info.chr),2,0);
    I_signif = sels{i}.GSS(tf) >= signif_thresh;
    signif(loc(I_signif),i) = 1;
    
end
save('combined_GSS_v2.mat', 'sels', 'signif');
%%
load('combined_GSS_v2.mat');
[tf,I] = ismemberCellRow(sels{1}.snp_id, sels{2}.snp_id);
sels{1} = cutStruc2Inds(sels{1}, find(tf));
%%
filein = 'SuppTable1_celUK1_5359_VeP_list_120219.csv';
supp = textscan(fopen(filein, 'r'), '%d %s %s %f %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 1);
supp  = {supp{1:3} [supp{4:8}]};
supp_sel = cell2struct(supp, {'rank', 'rs_1', 'rs_2', 'gss'}, 2);

supp_snp_id = [supp_sel.rs_1, supp_sel.rs_2];
val_snp_id = sels{1}.snp_id;

[supp_tf, val_I]=ismemberCellRow(supp_snp_id,val_snp_id) ;
val_I = val_I(supp_tf);

%%
cohort_i=1;
supp_I = find(supp_tf);
diff = supp_sel.gss(supp_tf,cohort_i)-sels{cohort_i}.GSS(val_I);
[sort_diff,I]=sort(diff);
plot(sort_diff)


















