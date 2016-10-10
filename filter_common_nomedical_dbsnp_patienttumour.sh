#!/bin/bash

. /etc/profile.d/modules.sh
module load julyin/tabix/0.2.6

VCF_LIB=/share/ClusterShare/software/contrib/julyin/vcflib/bin
DBSNP_DIR=/share/ClusterShare/biodata/contrib/julyin
DBSNP="common_no_known_medical_impact_20160502"


#REF_FILE=/share/ClusterShare/biodata/contrib/gi/gatk-resource-bundle/2.8/b37/human_g1k_v37.fasta
REF_FILE=/share/ClusterShare/biodata/contrib/gi/gatk-resource-bundle/2.5/hg19/ucsc.hg19.fasta

#SAMPLE_VCF=/home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-G89controlDNA-vs-G89HFMTissueDNA/myAnalysis/results/passed.somatic.snvs.vcf
SAMPLE_VCF=$1

$OUT_VCF=/home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-G89controlDNA-vs-G89HFMTissueDNA/myAnalysis/results/passed.somatic.snvs.filtcommon.vcf
OUT_VCF=$2

start_time=$(date +%s)
echo "[][][] Job started at "
date


<<-COM
sed -n -e '/^#/p' $SAMPLE_VCF > sample_header.txt
sed -n -e '/^chr1/p' $SAMPLE_VCF | head -10 > sample_chr1.txt
cat sample_header.txt sample_chr1.txt > sample_chr1.vcf
COMdd

echo "[][][] Creating dbsnp common vcf file with hg19 chromosome naming"
sed -n -e '/^#/p' $DBSNP_DIR/$DBSNP.vcf > $DBSNP_DIR/dbsnp_header.txt
sed -n -e '/^#/!p' $DBSNP_DIR/$DBSNP.vcf > $DBSNP_DIR/dbsnp_mutations.txt
sed -e 's/^/chr/' $DBSNP_DIR/dbsnp_mutations.txt > $DBSNP_DIR/dbsnp_mutations.chrmod.txt
cat $DBSNP_DIR/dbsnp_header.txt $DBSNP_DIR/dbsnp_mutations.chrmod.txt > $DBSNP_DIR/$DBSNP.chrmod.vcf

echo "[][][] Creating sample gz and tbi file"
bgzip -c -f $SAMPLE_VCF > $SAMPLE_VCF.gz
tabix -p vcf $SAMPLE_VCF.gz

echo "[][][] Creating dbsnp gz and tbi file"
bgzip -c -f $DBSNP_DIR/$DBSNP.chrmod.vcf > $DBSNP_DIR/$DBSNP.chrmod.vcf.gz
tabix -p vcf $DBSNP_DIR/$DBSNP.chrmod.vcf.gz
ddCOM

echo "[][][] Get complement "
/share/ClusterShare/software/contrib/julyin/bcftools/bcftools isec \
    --complement \
    --prefix /home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-G89controlDNA-vs-G89HFMTissueDNA/ \
    $SAMPLE_VCF.gz $DBSNP_DIR/$DBSNP.chrmod.vcf.gz
COM


echo "[][][] Get complement "
/share/ClusterShare/software/contrib/julyin/bcftools/bcftools isec \
    --prefix /home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-G89controlDNA-vs-G89HFMTissueDNA/passed.snvs.intersection \
    $SAMPLE_VCF.gz $DBSNP_DIR/$DBSNP.chrmod.vcf.gz



end_time=$(date +%s)
run_time=$(($end_time-$start_time))
((sec=run_time%60, run_time/=60, min=run_time%60, hrs=run_time/60))
timestamp=$(printf "[][][] Total time elapsed for script execution was: %d hours, %02d minutes, %02d seconds" $hrs $min $sec)
echo $timestamp


