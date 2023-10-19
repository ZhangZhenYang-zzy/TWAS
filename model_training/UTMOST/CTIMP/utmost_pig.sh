#!/bin/bash

# 需要把一个基因的多个表达量结果放到一个目录下
PLINK=/disk191/zzy/software/plink

SPECIE='pig'
TISSUE='Morula'
PRE_GENO=$TISSUE
GEX="/disk202/chenzt/.0_TWAS/Pig/expression_tmm_inv/${TISSUE}.expr_tmm_inv.bed.gz"
#gunzip -c $GEX > $TISSUE.expr_tmm_inv.bed
PRE_GEXP="$TISSUE.expr_tmm_inv.bed"

BATCH_START=0
BATCH_END=$(wc -l<$TISSUE.expr_tmm_inv.bed)
NR="${BATCH_START}_${BATCH_END}"
mkdir --parents $SPECIE/tmp
mkdir --parents $SPECIE/tmp/$TISSUE
mkdir --parents $SPECIE/tmp/$TISSUE/$NR


cat $PRE_GEXP | awk -vs=$BATCH_START -ve=$BATCH_END 'NR > s && NR <= e' |  while read PARAM; do

###取出基因左右1MB
CHR=`echo $PARAM | awk '{ print $1 }'`
P0=`echo $PARAM | awk '{ p=$2 - 1e6; if(p<0) p=0; print p; }'`
P1=`echo $PARAM | awk '{ print $3 + 1e6 }'`
GNAME=`echo $PARAM | awk '{ print $4 }'`
OUT="$SPECIE/tmp/$TISSUE/$NR/$PRE_GEXP.$GNAME"

echo $PARAM | tr ' ' '\n' | tail -n+5 | paste $PRE_GEXP.ID - > $OUT.pheno

$PLINK --bfile $PRE_GENO --pheno $OUT.pheno  --recode A --out $OUT --chr $CHR --from-bp $P0 --to-bp $P1



done







X=/disk191/zzy/.0_TWAS/2_software/CTIMP/example/X.txt
Y_folder=/disk191/zzy/.0_TWAS/2_software/CTIMP/example/Y_folder/
info=/disk191/zzy/.0_TWAS/2_software/CTIMP/example/info.txt
ntune=5 #number of grids for each tuning parameter
output_path=output/
mkdir ${output_path}
output_prefix=test # prefix of output files

Rscript /disk191/zzy/.0_TWAS/2_software/CTIMP/main.R ${X} ${info} ${Y_folder} ${ntune} ${output_prefix} ${output_path}
