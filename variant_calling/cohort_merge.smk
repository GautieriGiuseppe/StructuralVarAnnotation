# ==============================================================================
#  CONSOLIDATION
# ==============================================================================

rule octopusv_sample_merge_chm13:
    input:
        native_cutesv   = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_cuteSV.vcf',
        native_delly    = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_delly.vcf',
        native_sniffles = config['output'] + '/{batch}/{sample}/03.variant_calling/{sample}_chm13_sniffles.vcf',
        lifted_cutesv   = config['output'] + '/{batch}/{sample}/03.variant_calling/liftover/cuteSV_{sample}_grch38-to-chm13.vcf',
        lifted_delly    = config['output'] + '/{batch}/{sample}/03.variant_calling/liftover/delly_{sample}_grch38-to-chm13.vcf',
        lifted_sniffles = config['output'] + '/{batch}/{sample}/03.variant_calling/liftover/sniffles_{sample}_grch38-to-chm13.vcf'
    output:
        # The Union
        union_svcf = config['output'] + '/{batch}/{sample}/03.variant_calling/consolidated/{sample}_chm13_consolidated-union.svcf',
        union_vcf  = config['output'] + '/{batch}/{sample}/03.variant_calling/consolidated/{sample}_chm13_consolidated-union.vcf',
        # The Intersection & Plot 
        intersect_svcf = config['output'] + '/{batch}/{sample}/03.variant_calling/consolidated/{sample}_chm13_merged-intersect.svcf',
        upset_plot     = config['output'] + '/{batch}/{sample}/03.variant_calling/consolidated/{sample}_upset_plot.png'
    conda: "envs/octopus.yml"
    resources: mem_mb=config['hm'], time=config['ht']
    shell:
        '''
        TEMP_DIR="$(dirname {output.union_svcf})/temp_octopus_{wildcards.sample}"
        mkdir -p $TEMP_DIR

        # Standardize 6 files
        octopusv correct {input.native_cutesv} $TEMP_DIR/nc.svcf
        octopusv correct {input.native_delly} $TEMP_DIR/nd.svcf
        octopusv correct {input.native_sniffles} $TEMP_DIR/ns.svcf
        octopusv correct {input.lifted_cutesv} $TEMP_DIR/lc.svcf
        octopusv correct {input.lifted_delly} $TEMP_DIR/ld.svcf
        octopusv correct {input.lifted_sniffles} $TEMP_DIR/ls.svcf
        
        # Generate Union 
        octopusv merge -i $TEMP_DIR/*.svcf -o {output.union_svcf} --union
        octopusv svcf2vcf -i {output.union_svcf} -o {output.union_vcf}

        # Generate Intersection and Plots
        octopusv merge -i $TEMP_DIR/*.svcf \
            -o {output.intersect_svcf} --intersect \
            --upsetr --upsetr-output {output.upset_plot}

        rm -rf $TEMP_DIR
        '''

rule octopusv_cohort_merge_chm13:
    input:
        svcfs = [config['output'] + f"/{b}/{s}/03.variant_calling/consolidated/{s}_chm13_consolidated-union.svcf" for b, s in VALID_PAIRS]
    output:
        cohort_svcf = config['output'] + '/cohort_results/CHM13_final_cohort-union.svcf',
        cohort_vcf  = config['output'] + '/cohort_results/CHM13_final_cohort-union.vcf'
    conda: "envs/octopus.yml"
    resources: mem_mb=config['hm'], time=config['ht']
    shell:
        '''
        octopusv merge -i {input.svcfs} -o {output.cohort_svcf} --union
        octopusv svcf2vcf -i {output.cohort_svcf} -o {output.cohort_vcf}
        '''

# ==============================================================================
#  FORCE-CALLING (SNIFFLES2)
# ==============================================================================

rule sniffles2_genotype_chm13:
    """
    Genotypes each sample BAM against the consolidated cohort site list.
    """
    input:
        bam = config['output'] + '/{batch}/{sample}/01.align/chm13/{sample}.srt.bam',
        ref = config['reference']['chm13'],
        sites = config['output'] + '/cohort_results/CHM13_final_cohort-union.vcf'
    output:
        # We output to .vcf first as Sniffles2 native output
        vcf = config['output'] + '/{batch}/{sample}/04.force_calling/{sample}_chm13_genotyped.vcf'
    conda: "envs/sniffles.yml"
    threads: config['hc']
    resources:
        mem_mb=config['hm'],
        time=config['ht']
    shell:
        '''
        sniffles --input {input.bam} \
                 --reference {input.ref} \
                 --genotype-vcf {input.sites} \
                 --vcf {output.vcf} \
                 --sample {wildcards.sample} \
                 --threads {threads}
        '''

rule compress_and_index_genotypes:
    """
    Convert VCF to BCF for efficient merging fixing the missing Sniffles2 header tags
    """
    input:
        vcf = config['output'] + '/{batch}/{sample}/04.force_calling/{sample}_chm13_genotyped.vcf'
    output:
        bcf = config['output'] + '/{batch}/{sample}/04.force_calling/{sample}_chm13_genotyped.bcf',
        csi = config['output'] + '/{batch}/{sample}/04.force_calling/{sample}_chm13_genotyped.bcf.csi'
    conda: "envs/snakemake.yml"
    threads: 4
    resources: mem_mb=config['hm'], time=config['ht']
    shell:
        '''
        # Create header fix file
        echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">' > {wildcards.sample}_header_fix.txt
        echo '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of high-quality reference reads">' >> {wildcards.sample}_header_fix.txt
        echo '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of high-quality variant reads">' >> {wildcards.sample}_header_fix.txt

        # Step 1: Annotate ,sort and convert to BCF
        bcftools annotate -h {wildcards.sample}_header_fix.txt {input.vcf} | \
        bcftools sort -T {wildcards.sample}_sort_tmp -O b -o {output.bcf}

        # Step 2: Index the newly created BCF
        bcftools index --threads {threads} {output.bcf}

        # Clean up
        rm {wildcards.sample}_header_fix.txt
        '''

# ==============================================================================
#  COHORT MERGING
# ==============================================================================

rule merge_genotyped_cohort:
    """
    Merges all individual genotyped BCFs into a single population matrix.
    """
    input:
        bcfs = [config['output'] + f"/{b}/{s}/04.force_calling/{s}_chm13_genotyped.bcf" for b, s in VALID_PAIRS]
    output:
        vcf = config['output'] + '/NP057/cohort_results/CHM13_final_genotyped_matrix.vcf.gz'
    params:
        # Pass the sample IDs in the exact same order as the input BCFs
        sample_ids = [s for b, s in VALID_PAIRS]
    conda: "envs/snakemake.yml"
    resources: mem_mb=config['hm'], time=config['ht']
    shell:
        '''
        # 1. Create a temporary sample list file
        # This matches the order of the 'bcfs' input list exactly
        rm -f samples_to_rename.txt
        for s in {params.sample_ids}; do
            echo $s >> samples_to_rename.txt
        done

        # 2. Merge with --force-samples to bypass the 'SAMPLE' error
        # We pipe to reheader to apply the real names immediately
        bcftools merge --force-samples -m id -O u {input.bcfs} | \
        bcftools reheader -s samples_to_rename.txt | \
        bcftools view -O z -o {output.vcf}

        # 3. Index the final matrix
        bcftools index -t {output.vcf}

        # 4. Clean up
        rm samples_to_rename.txt
        '''