"""
Convert Delly sequence-based SV alleles to symbolic alleles.
Usage: python delly_to_symbolic.py input.vcf output.vcf
"""
import sys

in_path, out_path = sys.argv[1], sys.argv[2]

with open(in_path) as fin, open(out_path, 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue

        fields = line.strip().split('\t')
        if len(fields) < 8:
            fout.write(line)
            continue

        ref  = fields[3]
        alt  = fields[4]
        info = fields[7]
        pos  = int(fields[1])

        # Parse INFO into dict preserving flags (no-value entries)
        info_parts = info.split(';')
        info_dict  = {}
        info_flags = []
        for part in info_parts:
            if '=' in part:
                k, v = part.split('=', 1)
                info_dict[k] = v
            else:
                info_flags.append(part)

        svtype = info_dict.get('SVTYPE', '')
        end    = info_dict.get('END', '')

        # Convert sequence-based SVs to symbolic alleles
        # Only if: known SV type + long REF + ALT not already symbolic
        if svtype in ('DEL', 'DUP', 'INV') and len(ref) > 1 and not alt.startswith('<'):
            # Compute SVLEN from END - POS
            if end:
                svlen = int(end) - pos
            else:
                svlen = len(ref) - 1   # fallback: infer from REF length

            # DEL SVLEN is conventionally negative
            info_dict['SVLEN'] = str(-svlen) if svtype == 'DEL' else str(svlen)

            # Replace REF with anchor base, ALT with symbolic allele
            fields[3] = ref[0]
            fields[4] = f'<{svtype}>'

            # Rebuild INFO preserving original field order
            new_info = info_flags + [f'{k}={v}' for k, v in info_dict.items()]
            fields[7] = ';'.join(new_info)

        fout.write('\t'.join(fields) + '\n')