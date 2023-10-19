X=/disk191/zzy/.0_TWAS/2_software/CTIMP/example/X.txt
Y_folder=/disk191/zzy/.0_TWAS/2_software/CTIMP/example/Y_folder/
info=/disk191/zzy/.0_TWAS/2_software/CTIMP/example/info.txt
ntune=5 #number of grids for each tuning parameter
output_path=output/
mkdir ${output_path}
output_prefix=test # prefix of output files

Rscript /disk191/zzy/.0_TWAS/2_software/CTIMP/main.R ${X} ${info} ${Y_folder} ${ntune} ${output_prefix} ${output_path}