
# coding: utf-8

# Choose maf_file and desired_genes.
# 
# Currently variant types are considered potentially damaging. Filtering for chromosome type is currently commented out

# In[10]:

import re
import glob, os
from csv import DictReader, DictWriter


# Pick input and output file and list of desired genes.

# In[11]:

#maf_dir = "/home/julyin/genome_gbm_variant_calling/strelka_oncotator/strelka-oncotator-G53.normal.vs.G52.primary/"
maf_dir = "/home/julyin/genome_gbm_variant_calling/strelka_oncotator/strelka-oncotator-G53.normal.vs.G53.primary/"
os.chdir(maf_dir)


# with open('file_list.txt') as f:
#     files = f.read().splitlines()
#     
# print files

# In[12]:

with open('/home/julyin/analysis/multicentric_G52_G53/strelka/hg19_bam_pipeline_30x/gbm_gene_list.txt') as f:
    desired_genes = f.read().splitlines()
    
print desired_genes


# Desired variant damaging level and chromosomes.


variant_type_desired = dict.fromkeys(['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', \
        'In_Frame_Del', 'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site', \
        'Nonstop_Mutation', 'Start_Codon_Del', 'Start_Codon_Ins', 'Stop_Codon_Del', 'Stop_Codon_Ins', 'intergenic_variant'], 0)           
#chromosomes = range(1, 23) + ['X', 'Y']
#chromosome_list = set([str(x) for x in chromosomes])


# In[14]:

for mafFileName in glob.glob("headerless.*.maf"):
    
    print "Input: " + mafFileName
    outFileName = "WHOgenes." + ".".join(mafFileName.split(".")[1:])
    print  "Output: " + outFileName
    

    
    with open(mafFileName, 'rU') as inFile, \
        open(outFileName, "wb") as outFile:
        
        reader = DictReader(inFile, dialect='excel-tab')
        
        writer = DictWriter(f=outFile, fieldnames=reader.fieldnames, delimiter='\t')
        #writer.writerow(reader.fieldnames)
        writer.writeheader()
        
        for row in reader:   
               
            chr = row['Chromosome']
            start = row['Start_position']
            end = row['End_position']
            type = row['Variant_Type']
            gene = row['Hugo_Symbol']
            mutation_class = row['Variant_Classification']
                
            #if chr not in chromosome_list:          
            #    continue            
            if row['Matched_Norm_Sample_Barcode'] == 'NORMAL':          
                continue 
            
            if gene in desired_genes and mutation_class in variant_type_desired:
                #print row            
                writer.writerow(row)
                
    print "DONE CREATING OUTPUT: " + outFileName
    
print "DONE ALL!"





