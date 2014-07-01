function snp_info=readBIM(file_in, delim)
%Read in a BIM file
%Assume standard format of 6 columns

     bim_fields = {'chr', 'snp_id', 'pos', 'allele1', 'allele2'};
     bim = textscan(fopen(file_in,'r'), '%d %s %d %d %s %s', 'delimiter', delim);
     snp_info = cell2struct(bim(:,[1:2,4:6]), bim_fields, 2);     