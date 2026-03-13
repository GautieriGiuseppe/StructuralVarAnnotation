# ==============================================================================
#  VARIANT CALLING RULES
# ==============================================================================

# --- SNIFFLES ---
rule sniffles_grch38:
    input:  sample = config['output'] + '/{batch}/{sample}/01.align/grch38/{sample}.srt.bam'
    output: vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_grch38_sniffles.vcf',
            snf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_grch38_sniffles.snf'
    conda: "envs/sniffles.yml"
    threads: config['mc']
    resources: mem_mb=config['mm'], time=config['ht']
    shell: "sniffles --input {input.sample} --vcf {output.vcf} --snf {output.snf}"

rule sniffles_chm13:
    input:  sample = config['output'] + '/{batch}/{sample}/01.align/chm13/{sample}.srt.bam'
    output: vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_sniffles.vcf',
            snf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_sniffles.snf'
    conda: "envs/sniffles.yml"
    threads: config['mc']
    resources: mem_mb=config['mm'], time=config['ht']
    shell: "sniffles --input {input.sample} --vcf {output.vcf} --snf {output.snf}"

# --- DELLY ---
rule delly_grch38:
    input:  reference = config['reference']['grch38'],
            sample = config['output'] + '/{batch}/{sample}/01.align/grch38/{sample}.srt.bam'
    output: config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_grch38_delly.bcf'
    conda: "envs/delly.yml"
    threads: config['mc']
    resources: mem_mb=config['mm'], time=config['ht']
    shell: "./delly_static lr -y ont -g {input.reference} {input.sample} -o {output}"

rule delly_chm13:
    input:  reference = config['reference']['chm13'],
            sample = config['output'] + '/{batch}/{sample}/01.align/chm13/{sample}.srt.bam'
    output: config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_delly.bcf'
    conda: "envs/delly.yml"
    threads: config['lc']
    resources: mem_mb=config['mm'], time=config['ht']
    shell: "./delly_static lr -y ont -g {input.reference} {input.sample} -o {output}"

# --- CUTESV ---
rule cuteSV_grch38:
    input:  bam = config['output'] + '/{batch}/{sample}/01.align/grch38/{sample}.srt.bam',
            ref = config['reference']['grch38']
    output: vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_grch38_cuteSV.vcf'
    conda: "envs/cutesv.yml"
    threads: config['hc']
    resources: mem_mb=config['hm'], time=config['ht']
    shell:
        '''
        TEMP_DIR="$(dirname {output.vcf})/temp_cutesv_grch38_{wildcards.sample}"
        mkdir -p $TEMP_DIR
        cuteSV --threads {threads} --sample {wildcards.sample} --genotype {input.bam} {input.ref} {output.vcf} $TEMP_DIR
        rm -rf $TEMP_DIR
        '''

rule cuteSV_chm13:
    input:  
        bam = config['output'] + '/{batch}/{sample}/01.align/chm13/{sample}.srt.bam',
        ref = config['reference']['chm13'] 
    output: 
        vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_cuteSV.vcf'
    conda: "envs/cutesv.yml"
    threads: config['hc']
    resources: mem_mb=config['vhm'], time=config['vht']
    shell:
        '''
        TEMP_DIR="$(dirname {output.vcf})/temp_cutesv_chm13_{wildcards.sample}"
        rm -rf $TEMP_DIR 
        mkdir -p $TEMP_DIR
        
        cuteSV --threads {threads} --sample {wildcards.sample} --genotype {input.bam} {input.ref} {output.vcf} $TEMP_DIR
        
        rm -rf $TEMP_DIR
        '''

# --- BCFTOOLS ---
rule bcf_to_vcf:
    input:  "{path}/{sample}_{ref}_{tool}.bcf"
    output: "{path}/{sample}_{ref}_{tool}.vcf"
    wildcard_constraints:
        # This tells Snakemake: ONLY use this rule if the 'tool' is actually delly
        # This prevents it from trying to 'convert' liftover or cohort files incorrectly.
        tool = "delly"
    conda: "envs/snakemake.yml"
    resources: mem_mb=config['hm'], time=config['ht']
    shell: "bcftools view -O v {input} > {output}"

# ==============================================================================
#  LIFTOVER 
# ==============================================================================

rule crossmap_liftover:
    input:
        # We specify the vcf here. Because of our constraint above, 
        # Snakemake now knows it MUST run bcf_to_vcf first if the tool is delly.
        vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_{source}_{tool}.vcf',
        chain = config['chain_file'],
        target_ref = lambda wildcards: config['reference'][wildcards.target]
    output:
        # Changed the separator here to a dash '-' for the liftover part 
        # to prevent further underscore confusion
        vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/liftover/{tool}_{sample}_{source}-to-{target}.vcf'
    conda: "envs/crossmap.yml"
    threads: 8
    resources: mem_mb=config['hm'], time=config['ht']
    shell: "CrossMap vcf {input.chain} {input.vcf} {input.target_ref} {output.vcf}"