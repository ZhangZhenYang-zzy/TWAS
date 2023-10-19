# Training model for UTMOST

1_getGenelist.py: get the gene name list from all the bed files.

2_getExpressionEachtissue.py: extract the gene expression level from all the tissues.

3_ctimp.py: prepare the genotype and run CTIMP. run 3_ctimp.py start end. The start and end is the row number in the gene list (output by 1_getGenelist.py ), which will allow the users run the CTIMP in the different servers in the same time. 

4.get_UTMOST.corr.R: calculate the rsq for each gene and tissue.

5.create_table.R and 6.make_sqlite_db.py: come from https://github.com/Joker-Jerome/utmost_update

