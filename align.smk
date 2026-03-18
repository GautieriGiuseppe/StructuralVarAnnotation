import os

samples_file = config['samples']
print(f"Loading samples from: {samples_file}")

SAMPLES = []
BATCHES = set()
ubam_map = {}

with open(samples_file, 'r') as f:
    # Skip header
    header = f.readline().strip().split('\t')
    
    # Find column indices
    sample_idx = header.index('sample_id')
    ubam_idx = header.index('ubam')
    batch_idx = header.index('batch_id')
    
    # Read each line
    for line in f:
        if line.strip():  # Skip empty lines
            fields = line.strip().split('\t')
            sample = fields[sample_idx]
            ubam = fields[ubam_idx]
            batch = fields[batch_idx]
            
            SAMPLES.append(sample)
            BATCHES.add(batch)
            ubam_map[(batch, sample)] = ubam
            print(f"  {batch}/{sample} -> {ubam}")

SAMPLES = list(set(SAMPLES))  
BATCHES = list(BATCHES)

print(f"Loaded {len(SAMPLES)} samples: {SAMPLES}")
print(f"Batches: {BATCHES}")

# Define rule all
rule all:
    input:
        # GRCh38 alignments
        [f"{config['output']}/{batch}/{sample}/01.align/grch38/{sample}.srt.bam" 
         for batch in BATCHES for sample in SAMPLES],
        # CHM13 alignments
        [f"{config['output']}/{batch}/{sample}/01.align/chm13/{sample}.srt.bam" 
         for batch in BATCHES for sample in SAMPLES],

rule minimap2_GRCh38:
    input:
        reference=config['reference']['grch38'],
        sample=lambda wildcards: ubam_map.get((wildcards.batch, wildcards.sample), 
                                             f"{config['input']}/{wildcards.batch}/{wildcards.sample}.ubam")
    output:
        f"{config['output']}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam"
    threads:
        config['mc']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['hm'],
        time_min=lambda wildcards, attempt: attempt * config['vht']
    params:
        read_group=r"'@RG\tID:{sample}\tPL:ONT\tSM:{sample}'",
        prefix=f"{config['output']}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}"
    shell:
        '''
        samtools fastq \
            -T MM,ML \
            {input.sample} | \
        minimap2 \
            -a \
            -x map-ont \
            -t {threads} \
            --MD \
            -y \
            --rmq=yes \
            --cs \
            -R {params.read_group} \
            {input.reference} \
            - | \
        samtools sort \
            -m 3G \
            -@ {threads} \
            --write-index \
            -o {output} \
            -T {params.prefix} \
            -
        '''

rule minimap2_chm13:
    input:
        reference=config['reference']['chm13'],
        sample=lambda wildcards: ubam_map.get((wildcards.batch, wildcards.sample), 
                                             f"{config['input']}/{wildcards.batch}/{wildcards.sample}.ubam")
    output:
        f"{config['output']}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam"
    threads:
        config['mc']
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['hm'],
        time_min=lambda wildcards, attempt: attempt * config['vht']
    params:
        read_group=r"'@RG\tID:{sample}\tPL:ONT\tSM:{sample}'",
        prefix=f"{config['output']}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}"
    shell:
        '''
        samtools fastq \
            -T MM,ML \
            {input.sample} | \
        minimap2 \
            -a \
            -x map-ont \
            -t {threads} \
            --MD \
            -y \
            --rmq=yes \
            --cs \
            -R {params.read_group} \
            {input.reference} \
            - | \
        samtools sort \
            -m 3G \
            -@ {threads} \
            --write-index \
            -o {output} \
            -T {params.prefix} \
            -
        '''
