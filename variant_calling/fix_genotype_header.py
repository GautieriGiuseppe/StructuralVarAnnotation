"""
Fix Sniffles force-called VCF header:
  - Inherits 8-sample header from SURVIVOR cohort VCF
  - Has only 1 data column
  - Missing FORMAT tags (GQ, DR, DV)

Usage: python3 fix_genotype_header.py input.vcf output.vcf sample_name
"""
import sys

in_path, out_path, sample = sys.argv[1], sys.argv[2], sys.argv[3]

FORMAT_TAGS = [
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n',
    '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of high-quality reference reads">\n',
    '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of high-quality variant reads">\n',
]

with open(in_path) as fin, open(out_path, 'w') as fout:
    meta = []
    body = []
    for line in fin:
        if line.startswith('##'):
            meta.append(line)
        elif line.startswith('#CHROM'):
            pass  # skip old column header — will write new one
        else:
            body.append(line)

    # Add missing FORMAT tags
    existing = ''.join(meta)
    for tag in FORMAT_TAGS:
        tag_id = tag.split(',')[0]  # e.g. '##FORMAT=<ID=GQ'
        if tag_id not in existing:
            meta.append(tag)

    # Write fixed header
    fout.writelines(meta)
    fout.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', sample]) + '\n')
    fout.writelines(body)