gwis-stats
==========

A demonstration of the novel statistics implemented in 
the paper 
"Goudey (2013) - GWIS - model-free, fast and exhaustive search for epistatic interactions in case-control GWAS'



==How to run==


The main function is runGssAnalysis(plink_files, snp_list, GSS_type)

where 
plink_files is the filename of the PLINK binary files (.bed, .bim, .fam) without the final extension. This means all three files 
must share the same name. 

snp_list is a list of snp identifiers corresponding to those in the .bim file 
Can either be
 - A 1D list of identifiers. We will then exhaustively search over all possible pairs
    of these SNPs. 
 - A list of SNP pairs. Each of these will be explicityly evaluated. 
 
 
