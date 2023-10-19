#!/bin/bash

TISSUE=$1
BATCH_START=$2
BATCH_END=$3

#echo $TISSUE

GCTA="/path/to/gcta_nr_robust"
PLINK="/path/to/2_software/plink"
PLINK2="/path/to/plink2"
GEMMA="/path/to/gemma-0.98.1"


# LDREF="/disk201/chenzt/.0_TWAS/FUSION/Pig"
# TISSUE='Morula'
# GEX="/disk202/chenzt/.0_TWAS/Pig/expression_tmm_inv/${TISSUE}.expr_tmm_inv.bed.gz"
#gunzip -c $GEX > $TISSUE.expr_tmm_inv.bed
PRE_GEXP="$TISSUE.expr_tmm_inv.bed"

#cat $PRE_GEXP | head -n1 | tr '\t' '\n' | tail -n+5 | awk '{ print $1,$1 }' > ${PRE_GEXP}.ID

#$PLINK --vcf /disk201/chenzt/.0_TWAS/Pig/vcf_each_tissue/$TISSUE.filtered_maf0.05_mac6.geno.vcf.gz --double-id --make-bed --out /disk201/chenzt/.0_TWAS/FUSION/Pig/$TISSUE

PRE_GENO=$TISSUE
OUT_DIR="./WEIGHTS"

NR="${BATCH_START}_${BATCH_END}"
mkdir --parents tmp/$NR
mkdir --parents hsq/$NR
mkdir --parents out/$NR
mkdir $OUT_DIR
mkdir output
# ln -s /mnt/test/zzy/1_data_pig/tmp output

cat $PRE_GEXP | awk -vs=$BATCH_START -ve=$BATCH_END 'NR > s && NR <= e' |  while read PARAM; do
CHR=`echo $PARAM | awk '{ print $1 }'`


P0=`echo $PARAM | awk '{ print $3 - 1e6 }'`
if (( $P0 < 0 ))
then
	P0=0;
fi

P1=`echo $PARAM | awk '{ print $3 + 1e6 }'`
GNAME=`echo $PARAM | awk '{ print $4 }'`
OUT="tmp/$NR/$PRE_GEXP.$GNAME"

#  if exists the result
RESULTFILE="WEIGHTS/$TISSUE/$TISSUE.$GNAME.wgt.RDat"
if [ -f "$RESULTFILE" ]; then
    echo "$RESULTFILE exists."
    continue
fi


echo $PARAM | tr ' ' '\n' | tail -n+5 | paste $PRE_GEXP.ID - > $OUT.pheno
$PLINK --bfile $PRE_GENO --pheno $OUT.pheno --make-bed --out $OUT --keep $OUT.pheno --chr $CHR --from-bp $P0 --to-bp $P1
mkdir $OUT_DIR/$PRE_GENO
FINAL_OUT="$OUT_DIR/$PRE_GENO/$PRE_GENO.$GNAME"
Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --PATH_plink $PLINK2 --models top1,blup,bslmm,lasso --covar $TISSUE.covariates4Fusion.txt --crossval 5

cat $FINAL_OUT.hsq >> hsq/$NR.hsq
rm -f $FINAL_OUT.hsq $OUT.tmp.*
rm $OUT.*
done
