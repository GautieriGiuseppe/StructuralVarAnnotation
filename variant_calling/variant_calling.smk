import os

# ==============================================================================
#  PARSING CONFIG & METADATA
# ==============================================================================
samples_file = config['samples']
VALID_PAIRS = []

with open(samples_file, 'r') as f:
    header = f.readline().strip().split('\t')
    sample_idx = header.index('sample_id')
    batch_idx = header.index('batch_id')
    for line in f:
        if line.strip():
            fields = line.strip().split('\t')
            VALID_PAIRS.append((fields[batch_idx], fields[sample_idx]))

# ==============================================================================
#  INCLUDING THE TWO SMK FILES
# ==============================================================================
include: "sample_variant_calling.smk"
include: "cohort_merge.smk"

# ==============================================================================
#  RULE ALL
# ==============================================================================

rule all:
    input:
        # 1. Final Matrix
        config['output'] + '/NP057/cohort_results/CHM13_final_genotyped_matrix.vcf.gz',

        # 2. The Union 
        [config['output'] + f"/{b}/{s}/03.variant_calling/consolidated/{s}_chm13_consolidated-union.svcf" 
         for b, s in VALID_PAIRS],
        [config['output'] + f"/{b}/{s}/03.variant_calling/consolidated/{s}_chm13_consolidated-union.vcf" 
         for b, s in VALID_PAIRS],

        # 3. The Intersection 
        [config['output'] + f"/{b}/{s}/03.variant_calling/consolidated/{s}_chm13_merged-intersect.svcf" 
         for b, s in VALID_PAIRS],
        
        # 4. UpSet Plots
        [config['output'] + f"/{b}/{s}/03.variant_calling/consolidated/{s}_upset_plot.png" 
         for b, s in VALID_PAIRS]