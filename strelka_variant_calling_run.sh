#!/bin/bash
# Run Strelka on all GBM samples paired.

#qsub -q all.q -cwd -V -b y -j y -l h_vmem=8G -l mem_requested=8G -l pwbc=true -pe smp 8 -N strelka bash ./strelka_pdx_11SEP15.sh
#call by  ./script.sh [Work directory] [Path to normal] [Path to tumour]


start_time=$(date +%s)
echo "[][][] Job started at "
date


NORMAL_PATH=$1
TUMOUR_PATH=$2
NORMAL_tmp=$(basename $NORMAL_PATH)
TUMOUR_tmp=$(basename $TUMOUR_PATH)
#https://stackoverflow.com/questions/27658675/how-to-remove-last-n-characters-from-a-bash-variable-string
NORMAL=${NORMAL_tmp/.dedup.realn.bam/}
TUMOUR=${TUMOUR_tmp/.dedup.realn.bam/}

REF_FILE="/share/ClusterShare/biodata/contrib/gi/gatk-resource-bundle/2.8/hg19/ucsc.hg19.fasta"
#REF_FILE=$3

# example location where strelka was installed to:
#
STRELKA_INSTALL_DIR=/share/ClusterShare/software/contrib/julyin/strelka_workflow-1.0.15

echo "[][][] $NORMAL vs $TUMOUR strelka started at"
date

# example location where analysis will be run:
#
#WORK_DIR="/home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/strelka-$NORMAL-vs-$TUMOUR/"
WORK_DIR=$3
mkdir -p $WORK_DIR

# Step 1. Move to working directory:
#
cd $WORK_DIR

# Step 2. Copy configuration ini file from default template set to a
#         local copy, possibly edit settings in local copy of file:

cp $STRELKA_INSTALL_DIR/etc/strelka_config_bwa_default.ini strelka_config_bwa_default.ini

# Step 3. Configure:
#
$STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow.pl \
    --normal=$NORMAL_PATH \
    --tumor=$TUMOUR_PATH \
    --ref=$REF_FILE \
    --config=strelka_config_bwa_default.ini \
    --output-dir=./myAnalysis

# Step 4. Run Analysis
#         This example is run using 8 cores on the local host:
#
cd ./myAnalysis
make -j 8

echo "[][][] $NORMAL vs $TUMOUR strelka finished at"
date


end_time=$(date +%s)
run_time=$(($end_time-$start_time))
((sec=run_time%60, run_time/=60, min=run_time%60, hrs=run_time/60))
timestamp=$(printf "[][][] Total time elapsed for script execution was: %d hours, %02d minutes, %02d seconds" $hrs $min $sec)
echo $timestamp


