import os

# ==============================================================================
# 1. PARSING CONFIG & METADATA
# ==============================================================================

samples_file = config['samples']
print(f"Loading samples from: {samples_file}")

# We will store valid pairs here to drive 'rule all'
# Format: [(batch_id, sample_id), ...]
VALID_PAIRS = []

with open(samples_file, 'r') as f:
    # Skip header
    header = f.readline().strip().split('\t')
    
    # Find column indices
    sample_idx = header.index('sample_id')
    batch_idx = header.index('batch_id')
    
    # Read each line
    for line in f:
        if line.strip():  # Skip empty lines
            fields = line.strip().split('\t')
            sample = fields[sample_idx]
            batch = fields[batch_idx]
            
            # Store the valid pair
            VALID_PAIRS.append((batch, sample))

print(f"Loaded {len(VALID_PAIRS)} sample/batch pairs.")


# ==============================================================================
# 2. RULE ALL (The Driver)
# ==============================================================================

rule all:
    input:
        # We use list comprehensions to generate the exact targets needed
        # This prevents looking for Sample A in Batch B
        [f"{config['output']}/{batch}/{sample}/02.alignqc/grch38/{sample}.alfred.tsv.gz" for batch, sample in VALID_PAIRS],
        [f"{config['output']}/{batch}/{sample}/02.alignqc/grch38/{sample}.mosdepth.global.dist.txt" for batch, sample in VALID_PAIRS],
        [f"{config['output']}/{batch}/{sample}/02.alignqc/chm13/{sample}.alfred.tsv.gz" for batch, sample in VALID_PAIRS],
        [f"{config['output']}/{batch}/{sample}/02.alignqc/chm13/{sample}.mosdepth.global.dist.txt" for b, s in VALID_PAIRS]


# ==============================================================================
# 3. ANALYSIS RULES
# ==============================================================================

rule alfred_stats_grch38:
    input:
        reference=config['reference']['grch38'],
        sample=config['output'] + '/{batch}/{sample}/01.align/grch38/{sample}.srt.bam'
    output:
        config['output'] + '/{batch}/{sample}/02.alignqc/grch38/{sample}.alfred.tsv.gz'
    threads:
        config['lc']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['mm'],
        time=lambda wildcards, attempt: attempt * config['ht']
    params:
        sample='{sample}'
    shell:
        '''
        alfred qc \
            -r {input.reference} \
            -a {params.sample} \
            -o {output} \
            {input.sample}
        '''

rule mosdepth_stats_grch38:
    input:
        config['output'] + '/{batch}/{sample}/01.align/grch38/{sample}.srt.bam'
    output:
        config['output'] + '/{batch}/{sample}/02.alignqc/grch38/{sample}.mosdepth.global.dist.txt'
    threads:
        config['lc']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['mm'],
        time=lambda wildcards, attempt: attempt * config['mt']
    params:
        prefix=config['output'] + '/{batch}/{sample}/02.alignqc/grch38/{sample}'
    shell:
        '''
        mosdepth \
            {params.prefix} \
            {input} \
            -n \
            -x \
            --by 1000 
        '''

rule alfred_stats_chm13:
    input:
        reference=config['reference']['chm13'],
        sample=config['output'] + '/{batch}/{sample}/01.align/chm13/{sample}.srt.bam'
    output:
        config['output'] + '/{batch}/{sample}/02.alignqc/chm13/{sample}.alfred.tsv.gz'
    threads:
        config['lc']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['mm'],
        time=lambda wildcards, attempt: attempt * config['ht']
    params:
        sample='{sample}'
    shell:
        '''
        alfred qc \
            -r {input.reference} \
            -a {params.sample} \
            -o {output} \
            {input.sample}
        '''

rule mosdepth_stats_chm13:
    input:
        config['output'] + '/{batch}/{sample}/01.align/chm13/{sample}.srt.bam'
    output:
        config['output'] + '/{batch}/{sample}/02.alignqc/chm13/{sample}.mosdepth.global.dist.txt'
    threads:
        config['lc']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['mm'],
        time=lambda wildcards, attempt: attempt * config['mt']
    params:
        prefix=config['output'] + '/{batch}/{sample}/02.alignqc/chm13/{sample}'
    shell:
        '''
        mosdepth \
            {params.prefix} \
            {input} \
            -n \
            -x \
            --by 1000 
        '''