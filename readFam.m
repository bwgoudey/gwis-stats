function sample_info=readFam(file_in, delim)
%Read in a FAM file
%Assume standard format of 6 columns

     fam_fields = {'sample_id', 'indiv_id', 'pat_id', 'mat_id', 'sex', 'label'};
     bim = textscan(fopen(file_in,'r'), '%s %s %s %s %d %d', 'delimiter', delim);
     keep_fields = [1,2,5,6];
     
     sample_info = cell2struct(bim(:,keep_fields), fam_fields(keep_fields), 2);     