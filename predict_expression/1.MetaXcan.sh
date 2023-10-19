tissue=Muscle
vcf_genotype=using.vcf.gz
db=${tissue}.db
nohup python3 MetaXcan/software/Predict.py \
--model_db_path $db \
--vcf_genotypes $vcf_genotype \
--vcf_mode genotyped \
--prediction_output result/${tissue}_predict.txt \
--prediction_summary_output result_summary/${tissue}_summary.txt \
--verbosity 9 \
--throw