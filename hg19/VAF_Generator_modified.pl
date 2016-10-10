#!/usr/bin/perl -w

use strict; 
use warnings;

use FindBin qw($Bin); # Directory the script is in
use Sort::Naturally; # Allows natural sorting with nsort: http://search.cpan.org/~bingos/Sort-Naturally-1.03/lib/Sort/Naturally.pm

# This script takes VEP-annotated Strelka outputs for any number of Normal-Tumour pairs and outputs a CSV with each variant that is present in one or more of these comparisons. Information output for each variant includes the depth, quality, allele frequency and variant impact.

####################################



#my $normal_sample=$ARGV[0];
#my $tumour_sample=$ARGV[1];

my $first_vcf=$ARGV[0];
my $second_vcf=$ARGV[1];
my $user_out=$ARGV[2];

#print "Hello, $normal_sample $tumour_sample\n";

my @vcf_files;
#push (@vcf_files,$Bin."/bwa_strelka/41414_4_vs_414141A/passed.somatic.snvs.VEP.vcf");
#push (@vcf_files,$Bin."/bwa_strelka/41414_4_vs_41414_5/passed.somatic.snvs.VEP.vcf");
#push (@vcf_files,$Bin."/bwa_strelka/41414_4_vs_41414_2/passed.somatic.snvs.VEP.vcf");
#push (@vcf_files, "/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/G89controlDNA.vs.G89PDX.commonfilt.passed.snvs.vep.vcf");
#push (@vcf_files, "/home/julyin/genome_gbm_variant_calling/strelka_VEP_annotated/G89controlDNA.vs.G89HFMTissue.commonfilt.passed.snvs.vep.vcf");
push (@vcf_files, $first_vcf);
push (@vcf_files, $second_vcf);


#my $outputfile = "/home/julyin/analysis/multicentric_G52_G53/strelka/TEST.csv";
my $outputfile = $user_out;

#my $outputfile = "/home/julyin/genome_gbm_variant_calling/strelka/strelka-G53controlDNA-vs-G53TissueDNA/myAnalysis/results/vep_annotated_results/VAFs-G53controlDNA-vs-G53TissueDNA-passed.somatic.snvs.csv";

my %vcf_hash;
my %vaf;
my %unique_variants;
my %variant_counts;

####################################

# Hash the query vcfs with the VAFs for each variant

print "Hashing input VCFs.\n";

foreach my $vcf_file (@vcf_files) {
	hash_vcf($vcf_file);
}

print "Done.\n";

####################################

open(OUTPUT, ">$outputfile") or die "Cannot open file \"$outputfile\" for output.\n\n";

# Print a key to what each heading line is
for (my $i=0; $i<scalar(@vcf_files);$i++) {
	print OUTPUT "##".($i+1)."=".$vcf_files[$i]."\n";
}

print OUTPUT "Chromosome,Position,Effect,Gene(s),";
for (my $i=0; $i<scalar(@vcf_files);$i++) {
	print OUTPUT "SGT".($i+1).",QSS".($i+1).",VAF".($i+1)."_T1,VAF".($i+1)."_T2,DP".($i+1)."_T1,DP".($i+1)."_T2";
	if (($i+1) != scalar(@vcf_files)) { print OUTPUT ","; } # Only print a comma at the end if there are more VCF files to put in
}
print OUTPUT "\n";

foreach my $chr (nsort keys %unique_variants) {
	foreach my $pos (nsort keys $unique_variants{$chr}) {
		$variant_counts{$unique_variants{$chr}{$pos}}++; # Count the number of variants present in 1, 2, 3, etc VCF files for reporting later
		
		my $vep_annotation_flag = 0; # Create a flag for determining whether the vep annotation information has been printed already

		my $output = $chr.",".$pos.","; # Create the output line and populate it with the first 2 columns
		
		# Determine variant effects and genes affected using the VEP annotation and print this as the next 2 columns
		for (my $i=0; $i<scalar(@vcf_files);$i++) {
			if (defined $vcf_hash{$vcf_files[$i]}{$chr}{$pos}) {
				if ($vep_annotation_flag == 0) {
					$output .= vep_annotation($vcf_hash{$vcf_files[$i]}{$chr}{$pos}).",";
					
					$vep_annotation_flag = 1; # Set flag to true so further annotations are not printed
				}
			}
		}
		
		# Print the sets of columns for each variant file with sgt, qss, vaf and depth information
		for (my $i=0; $i<scalar(@vcf_files);$i++) {
			if (defined $vcf_hash{$vcf_files[$i]}{$chr}{$pos}) {
				$output .= sgt_qss($vcf_hash{$vcf_files[$i]}{$chr}{$pos}).",".vaf_depth($vcf_hash{$vcf_files[$i]}{$chr}{$pos}).",";

			} else {
				$output .= "N/A,0,0,0,0,0,";
			}	
		}
		
		chop($output); # Remove that last comma at the end of the output
		
		$output .= "\n";
		
		print OUTPUT $output;
	}
}

print "Number of variants present in X VCFs:\n";
foreach my $var_freq (nsort keys %variant_counts) {
	print $var_freq.": ".$variant_counts{$var_freq}."\n";
}

close OUTPUT;

####################################

exit;

####################################

sub hash_vcf {
	my ($query_vcf) = @_;

	-e $query_vcf or die "File \"$query_vcf\" does not exist.\n";
	
	open(QUERYVCF, $query_vcf) or die "Cannot open file \"$query_vcf\"\n\n"; 

	foreach my $line (<QUERYVCF>) {
		if ($line =~ /^([\.\w]*)\t([0-9]*)\t[\w\.]\t(\w)\t(\w)\t[\w\.]\t[\w]*\t(.*?)\t[\w:]*\t[0-9:,]*\t([0-9:,]*)/) { # Match every variant line in the VCF
			my $chr = $1;
			my $pos = $2;
			my $ref = $3;
			my $alt = $4;
			my $info = $5;
			my $tumor_depth_af = $6;
			
			$vcf_hash{$query_vcf}{$chr}{$pos} = "ALT=".$alt."\t".$info."\t".$tumor_depth_af; # Save variant info for calculations later
			
			if (!defined $unique_variants{$chr}{$pos}) { # If the current position has not been seen before, set its frequency to 1
				$unique_variants{$chr}{$pos} = 1;
			} else { # Otherwise increment the frequency by 1
				$unique_variants{$chr}{$pos}++;
			}
		} elsif ($line =~ /^#.*/) { # Match and ignore every info line in the VCF
			next;
		} else { # Catch all lines not matched by the regexes and print warning
			###### CURRENTLY IGNORING MULTI-ALLELIC SITES
			#print "Warning! Found a line in the query VCF file that can't be parsed with the regex:\n".$line."\n";
		}
	}
	close QUERYVCF;
}

####################################

sub sgt_qss { # Snipes the SGT and QSS fields only 
	my ($depth_af_block) = @_;
	my $genotype;
	my $qss;
	
	if ($depth_af_block =~ /SGT=([\w\-\>]*);/) {
		$genotype = $1;
	}
	
	if ($depth_af_block =~ /QSS=([0-9]*);/) {
		$qss = $1;
	}
	
	return $genotype.",".$qss;
}

####################################

sub vaf_depth { # Only processes the last block of the VCF holding AF and depth info
	my ($depth_af_block) = @_;
	
	my $alt_allele;
	
	if ($depth_af_block =~ /ALT=([A-Z])/) {
		$alt_allele = $1;
	} else {
		print "Can't find alt allele for the variant containing ".$depth_af_block."\n";
		exit;
	}
	
	$depth_af_block =~ s/^.*\t(.*)/$1/; # Only keep the last block holding the depth and AF of the tumour genome
	
	my @split_block = split(/:/, $depth_af_block);
	
	my @AU = split(/,/, $split_block[4]);
	my @CU = split(/,/, $split_block[5]);
	my @GU = split(/,/, $split_block[6]);
	my @TU = split(/,/, $split_block[7]);
	
	my $tier1_depth = $AU[0] + $CU[0] + $GU[0] + $TU[0];
	my $tier2_depth = $AU[1] + $CU[1] + $GU[1] + $TU[1];
		
	if ($alt_allele eq "A") {
		return $AU[0] / $tier1_depth.",".$AU[1] / $tier2_depth.",".$tier1_depth.",".$tier2_depth;
	} elsif ($alt_allele eq "C") {
		return $CU[0] / $tier1_depth.",".$CU[1] / $tier2_depth.",".$tier1_depth.",".$tier2_depth;
	} elsif ($alt_allele eq "T") {
		return $TU[0] / $tier1_depth.",".$TU[1] / $tier2_depth.",".$tier1_depth.",".$tier2_depth;
	} elsif ($alt_allele eq "G") {
		return $GU[0] / $tier1_depth.",".$GU[1] / $tier2_depth.",".$tier1_depth.",".$tier2_depth;
	} else {
		print "Warning! An alt allele was passed for VAF calculation that is not A, C, T or G. It was: ".$alt_allele."\n";
	}
}

####################################

sub vep_annotation { # Only processes the last block of the VCF holding AF and depth info
	my ($vep_annotation) = @_;
	
	my %impacts; # Stores impact category counts
	$impacts{"high"} = 0; $impacts{"moderate"} = 0; $impacts{"modifier"} = 0; $impacts{"low"} = 0; # Set impacts to 0 to avoid errors when printing later
	my @genes_affected; # Stores genes affected

	$vep_annotation =~ s/.*CSQ=(.*?)\t.*/$1/; # Only keep the first block holding the VEP annotations

	while ($vep_annotation =~ /(.*?)\|.*?\|.*?\|.*?\|(.*?)\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?\|.*?,*/g) {
		my $impacts = $1;
		my $gene = $2;
		
		if (!grep(/^$gene$/, @genes_affected)) { # If the gene is not already on the list of affected genes, add it
  			push (@genes_affected, $gene);
		}
		
		my @multi_impacts = split ("&", $impacts); 
		
		for (my $i=0; $i<scalar(@multi_impacts);$i++) {
			# HIGH
			if (grep(/^frameshift_variant$/, @multi_impacts) || grep(/^splice_acceptor_variant$/, @multi_impacts) || grep(/^splice_donor_variant$/, @multi_impacts) || grep(/^stop_lost$/, @multi_impacts) || grep(/^stop_gained$/, @multi_impacts) || grep(/^transcript_ablation$/, @multi_impacts) || grep(/^transcript_amplification$/, @multi_impacts) || grep(/^incomplete_terminal_codon_variant$/, @multi_impacts) || grep(/^NMD_transcript_variant$/, @multi_impacts) || grep(/^feature_elongation$/, @multi_impacts) || grep(/^feature_truncation$/, @multi_impacts)) {
				$impacts{"high"}++;
			# MODERATE
			} elsif (grep(/^coding_sequence_variant$/, @multi_impacts) || grep(/^inframe_insertion$/, @multi_impacts) || grep(/^inframe_deletion$/, @multi_impacts) || grep(/^missense_variant$/, @multi_impacts) || grep(/^splice_region_variant$/, @multi_impacts) || grep(/^mature_miRNA_variant$/, @multi_impacts)) {
				$impacts{"moderate"}++;
			# MODIFIER
			} elsif (grep(/^downstream_gene_variant$/, @multi_impacts) || grep(/^exon_variant$/, @multi_impacts) || grep(/^non_coding_exon_variant$/, @multi_impacts) || grep(/^upstream_gene_variant$/, @multi_impacts) || grep(/^downstream_gene_variant$/, @multi_impacts) || grep(/^intron_variant$/, @multi_impacts) || grep(/^regulatory_region_variant$/, @multi_impacts) || grep(/^3_prime_UTR_variant$/, @multi_impacts) || grep(/^5_prime_UTR_variant$/, @multi_impacts) || grep(/^non_coding_transcript_exon_variant$/, @multi_impacts) || grep(/^non_coding_transcript_variant$/, @multi_impacts) || grep(/^nc_transcript_variant$/, @multi_impacts)) {
				$impacts{"modifier"}++;
			# LOW
			} elsif (grep(/^initiator_codon_variant$/, @multi_impacts) || grep(/^stop_retained_variant$/, @multi_impacts) || grep(/^synonymous_variant$/, @multi_impacts) || grep(/^TFBS_ablation$/, @multi_impacts) || grep(/^TFBS_amplification$/, @multi_impacts) || grep(/^TF_binding_site_variant$/, @multi_impacts) || grep(/^regulatory_region_ablation$/, @multi_impacts) || grep(/^regulatory_region_amplification$/, @multi_impacts) || grep(/^intergenic_variant$/, @multi_impacts)) {
				$impacts{"low"}++;
			} else {
				print "WARNING: could not classify impact of variant with impact: ".$multi_impacts[$i]."\n";
			}
		}
	}
	
	my $genes_affected_list = join(";", @genes_affected);
		
	return "high:".$impacts{"high"}.";moderate:".$impacts{"moderate"}.";modifier:".$impacts{"modifier"}.";low:".$impacts{"low"}.",".$genes_affected_list;
}
