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
    resources:
        mem_mb = config['mm'],
        time   = config['ht']
    shell: "sniffles --input {input.sample} --vcf {output.vcf} --snf {output.snf} --threads {threads}"

rule sniffles_chm13:
    input:  sample = config['output'] + '/{batch}/{sample}/01.align/chm13/{sample}.srt.bam'
    output: vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_sniffles.vcf',
            snf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_sniffles.snf'
    conda: "envs/sniffles.yml"
    threads: config['mc']
    resources:
        mem_mb = config['mm'],
        time   = config['ht']
    shell: "sniffles --input {input.sample} --vcf {output.vcf} --snf {output.snf} --threads {threads}"

# --- DELLY ---
rule delly_grch38:
    input:  reference = config['reference']['grch38'],
            sample = config['output'] + '/{batch}/{sample}/01.align/grch38/{sample}.srt.bam'
    output: config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_grch38_delly.bcf'
    conda: "envs/delly.yml"
    threads: config['mc']
    resources:
        mem_mb = config['mm'],
        time   = config['ht']
    shell: "delly lr -y ont -g {input.reference} {input.sample} -o {output}"

rule delly_chm13:
    input:  reference = config['reference']['chm13'],
            sample = config['output'] + '/{batch}/{sample}/01.align/chm13/{sample}.srt.bam'
    output: config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_delly.bcf'
    conda: "envs/delly.yml"
    threads: config['lc']
    resources:
        mem_mb = config['mm'],
        time   = config['ht']
    shell: "delly lr -y ont -g {input.reference} {input.sample} -o {output}"

# --- CUTESV ---
rule cuteSV_grch38:
    input:  bam = config['output'] + '/{batch}/{sample}/01.align/grch38/{sample}.srt.bam',
            ref = config['reference']['grch38']
    output: vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_grch38_cuteSV.vcf'
    conda: "envs/cutesv.yml"
    threads: config['hc']
    resources:
        mem_mb = config['hm'],
        time   = config['ht']
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
    resources:
        mem_mb = config['vhm'],
        time   = config['vht']
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
        tool = "delly"
    conda: "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb = config['hm'],
        time   = config['ht']
    shell: "bcftools view -O v {input} > {output}"

# ==============================================================================
#  LIFTOVER
# ==============================================================================

rule delly_to_symbolic:
    input:
        vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_{ref}_delly.vcf'
    output:
        vcf = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_{ref}_delly_symbolic.vcf'
    conda: "envs/snakemake.yml"
    threads: 1
    resources:
        mem_mb = config['mm'],
        time   = config['ht']
    shell: "python3 delly_to_symbolic.py {input.vcf} {output.vcf}"
