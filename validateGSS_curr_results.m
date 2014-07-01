addpath(genpath('/home/bwgoudey/Research/gwis-pp/post_process'))
addpath(genpath('/home/bwgoudey/Research/gwis-pp/source'))
addpath(genpath('/home/bwgoudey/Research/gwis-pp/binaries'))

cohorts={'UK1b', 'UK2b', 'FIN', 'NL', 'IT'};
proj_dir='/media/ntfs_2TB/Research/gwis-pp/projects/cel/';
getPlinkName = @(cohort) sprintf('%s/data/cel%s_noX_HH/A_raw/cel%s_noX_HH', proj_dir, cohort,cohort);

getBadGenoName = @(cohort) sprintf('%s/data/cel%s_noX_HH/B_annt/cel%s_noX_HH_bad_geno.mat', proj_dir, cohort,cohort);

%% Load data
plink_files = getPlinkName(cohorts{1});
snp_info=readBim([plink_files '.bim'], '\t');

%%
filein = 'SuppTable1_celUK1_5359_VeP_list_120219.csv';
supp = textscan(fopen(filein, 'r'), '%d %s %s %f %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 1);
supp  = {supp{1:3} [supp{4:8}]};
supp_sel = cell2struct(supp, {'rank', 'rs_1', 'rs_2', 'gss'}, 2);
supp_sel.snp_id = [supp_sel.rs_1, supp_sel.rs_2];

%% Get signif in other cohorts. 
signif = zeros(length(supp_sel.snp_id),length(cohorts));

for i = 1:length(cohorts)
    plink_files = getPlinkName(cohorts{i});
    snp_info=readBim([plink_files '.bim'], '\t');
    
    sels{i} = runGssAnalysis(plink_files, supp_sel.snp_id, 'GSS');
    
    
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
cohort_i=1;
supp_I = find(supp_tf);
diff = supp_sel.gss(supp_tf,cohort_i)-sels{cohort_i}.GSS(val_I);
[sort_diff,I]=sort(diff);
plot(sort_diff)


















