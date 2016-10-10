#!/bin/bash


SCRIPT="/home/julyin/analysis/template_scripts/vep_annotation.sh"

<<-COM
INPUT="/home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-G53.normal.vs.G52.primary/myAnalysis/results/passed.somatic.snvs.vcf"
OUT="/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/strelka-vep-G53.normal.vs.G52.primary/vep.passed.somatic.snvs.vcf"
STATS="/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/strelka-vep-G53.normal.vs.G52.primary/vep.stats.html"
NAME="log_G52_vep_hg19"

bash $SCRIPT $INPUT $OUT $STATS 2>&1 | tee log_vep_annotation_G52.log
COM

INPUT="/home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-G53.normal.vs.G53.primary/myAnalysis/results/passed.somatic.snvs.vcf"
OUT="/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/strelka-vep-G53.normal.vs.G53.primary/vep.passed.somatic.snvs.vcf"
STATS="/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/strelka-vep-G53.normal.vs.G53.primary/vep.stats.html"
NAME="log_G53_vep_hg19"

bash $SCRIPT $INPUT $OUT $STATS 2>&1 | tee log_vep_annotation_G53.log


