# Run strelka on multicentric GBM
# 27 September 2016
# Julia X.M. Yin


configfile: "/home/julyin/analysis/template_scripts/config.yaml"

# Adjust these variables
SAMPLE_LIST = ["G52.primary", "G53.primary"]
NORMAL = "G53.control"
TUMOUR = "primary"
BAM_DIR = "/share/ScratchGeneral/julyin/wgs_hg19_bam_files/run/phase2/"
VCF_DIR = "/home/julyin/genome_gbm_variant_calling/multicentric_G52_G53/hg19_30x/strelka-snvs-indels/"
ONCO_DIR = "/home/julyin/genome_gbm_variant_calling/multicentric_G52_G53/hg19_30x/strelka_oncotator/"
VEP_DIR = "/home/julyin/genome_gbm_variant_calling/multicentric_G52_G53/hg19_30x/strelka_VEP_annotated/"
VAR_TYPE = ["snvs", "indels"]

rule all:
    input: 
        expand(ONCO_DIR + "strelka-oncotator-" + NORMAL + ".vs.{sample}/{sample}.{type}.oncotator.maf", sample=SAMPLE_LIST, type=VAR_TYPE)
        #"/home/julyin/genome_gbm_variant_calling/vaf/vaf.multicentric.G52.vs.G53.csv"


rule strelka:
    input:
        normal_bam = BAM_DIR + NORMAL + ".dedup.realn.bam",
        tumour_bam = BAM_DIR + "{sample}.dedup.realn.bam"
    output:
        snvs = VCF_DIR + "strelka-" + NORMAL + ".vs.{sample}/myAnalysis/results/passed.somatic.snvs.vcf",
        indels = VCF_DIR + "strelka-" + NORMAL + ".vs.{sample}/myAnalysis/results/passed.somatic.indels.vcf"
    params:
        outdir = VCF_DIR + "strelka-" + NORMAL + ".vs.{sample}"
    message:
        """ [][][] Strelka variant calling [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        bash /home/julyin/analysis/template_scripts/strelka_snakemake_template/strelka_variant_calling_run.sh {input.normal_bam} {input.tumour_bam} {params.outdir}
        """

rule oncotator:
    input:
        VCF_DIR + "strelka-" + NORMAL + ".vs.{sample}/myAnalysis/results/passed.somatic.{type}.vcf"
    output:
        ONCO_DIR + "strelka-oncotator-" + NORMAL + ".vs.{sample}/{sample}.{type}.oncotator.maf"
    message:
        """ [][][] Oncotator annotation [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
    source /share/ClusterShare/software/contrib/julyin/oncotator-venv/bin/activate
    {config[software][oncotator]} -i VCF \
            --db-dir={config[refs][onc_db]} \
            --canonical-tx-file={config[refs][onc_canonical_tx_file]} \
            --output_format=TCGAMAF \
            {input} {output} hg19
        """



rule vep_annotation:
    input:
        VCF_DIR + "strelka-" + NORMAL + ".vs.{sample}/myAnalysis/results/passed.somatic.{type}.vcf"
    output:
        vep = VEP_DIR + "strelka-vep-" + NORMAL + ".vs.{sample}/{sample}.{type}.vep.vcf",
        stats = VEP_DIR + "strelka-vep-" + NORMAL + ".vs.{sample}/{sample}.{type}.vep.stats.html"
    message:
        """ [][][] VEP annotation [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: 
        """
        module load marcow/perl/5.14.2
        module load vep/76
        module load samtools/0.1.19
        perl /share/ClusterShare/software/contrib/gi/vep/76/variant_effect_predictor.pl \
             --cache \
             --dir /share/ClusterShare/biodata/contrib/gi/vep \
             --port 3337 \
             --offline \
             --input_file {input} \
             --format vcf \
             --output_file {output.vep} \
             --vcf \
             --stats_file {output.stats} \
             --stats_text \
             --force_overwrite \
             --canonical \
             --fork 32 \
             --sift b \
             --polyphen b \
             --symbol \
             --numbers \
             --terms so \
             --biotype \
             --total_length \
             --plugin LoF,human_ancestor_fa:/share/ClusterShare/biodata/contrib/gi/LOFTEE/1.0/human_ancestor.fa.rz --fields Consequenc
        e,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,Feature_type,cDNA_position,CDS_
        position,Existing_variation,DISTANCE,STRAND,CLIN_SIG,LoF_flags,LoF_filter,LoF,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,
        CADD_raw,CADD_phred,Reliability_index
        """


rule vaf_extraction:
    input:
        expand(VEP_DIR + "strelka-vep-" + NORMAL + ".vs.{sample}/{sample}.snvs.vep.vcf", sample = SAMPLE_LIST)
        #g53 = VCF_DIR + "strelka-" + NORMAL + ".vs.G53.primary/myAnalysis/results/passed.somatic.snvs.vcf"
    output:
        "/home/julyin/genome_gbm_variant_calling/vaf/vaf.multicentric.G52.vs.G53.csv"
    message:
        """ [][][] VAF extraction [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        /opt/perl/bin/perl /home/julyin/analysis/multicentric_G52_G53/strelka/VAF_Generator_modified.pl \
            {input} {output}
            """


