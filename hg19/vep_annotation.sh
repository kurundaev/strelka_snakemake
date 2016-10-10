#!/bin/bash


SCRIPT="/home/julyin/analysis/template_scripts/vep_annotation.sh"


INPUT="/home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-G53.normal.vs.G52.primary/myAnalysis/results/passed.somatic.snvs.vcf"
OUT="/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/strelka-vep-G53.normal.vs.G52.primary/vep.passed.somatic.snvs.vcf"
STATS="/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/strelka-vep-G53.normal.vs.G52.primary/vep.stats.html"
NAME="log_G52_vep_hg19"

qsub -q all.q -cwd -b y -j y -l h_vmem=8G -l mem_requested=8G -N $NAME bash $SCRIPT $INPUT $OUT $STATS


INPUT="/home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-G53.normal.vs.G53.primary/myAnalysis/results/passed.somatic.snvs.vcf"
OUT="/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/strelka-vep-G53.normal.vs.G53.primary/vep.passed.somatic.snvs.vcf"
STATS="/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/strelka-vep-G53.normal.vs.G53.primary/vep.stats.html"
NAME="log_G53_vep_hg19"

qsub -q all.q -cwd -b y -j y -l h_vmem=8G -l mem_requested=8G -N $NAME bash $SCRIPT $INPUT $OUT $STATS


