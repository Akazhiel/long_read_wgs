"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
"""

import re


def reformat_nanomonsv(inp, out):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('##') and 'ID=TR' in line:
            filtered_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            new_DR = line.replace(
                'ID=TR,Number=1,Type=Integer,Description="The number of reads around the breakpoints"',
                'ID=DR,Number=1,Type=Integer,Description="# of reads supporting the reference allele."',
            )
            filtered_vcf.write(new_DR)
        elif line.startswith('##') and 'ID=VR' in line:
            new_DV = line.replace(
                'ID=VR,Number=1,Type=Integer,Description="The number of variant supporting reads determined in the validation realignment step"',
                'ID=DV,Number=1,Type=Integer,Description="# of reads supporting the variant allele."',
            )
            filtered_vcf.write(new_DV)
        elif line.startswith('#CHROM'):
            headers = (
                line.strip()
                .replace('TUMOR', 'TUMOR_nanomon')
                .replace('CONTROL', 'NORMAL_nanomon')
                .split('\t')
            )
            filtered_vcf.write(line.replace('TUMOR', 'TUMOR_nanomon').replace('CONTROL', 'NORMAL_nanomon'))
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if columns[headers.index('FILTER')] == 'PASS':
                Format = columns[headers.index('FORMAT')].replace('TR', 'DR').replace('VR', 'DV')
                Format = 'GT:' + Format
                Normal_Format = columns[headers.index('NORMAL_nanomon')].split(':')
                Normal_Format[0] = str(int(Normal_Format[0]) - int(Normal_Format[1]))
                Tumor_Format = columns[headers.index('TUMOR_nanomon')].split(':')
                Tumor_Format[0] = str(int(Tumor_Format[0]) - int(Tumor_Format[1]))
                Normal = './.:' + ':'.join(Normal_Format)
                Tumor = './.:' + ':'.join(Tumor_Format)
                filtered_vcf.write(
                    '{}\t{}\t{}\t{}\n'.format('\t'.join(columns[0:8]), Format, Tumor, Normal)
                )
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()


def reformat_svim(inp, out, columnid, qual):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('##') and 'ID=DP' in line:
            new_DR = line.replace(
                'ID=DP,Number=1,Type=Integer,Description="Read depth"',
                'ID=DR,Number=1,Type=Integer,Description="# reads supporting the reference allele."',
            )
            filtered_vcf.write(new_DR)
        elif line.startswith('##') and 'ID=AD' in line:
            new_DV = line.replace(
                'ID=AD,Number=R,Type=Integer,Description="Read depth for each allele"',
                'ID=DV,Number=R,Type=Integer,Description="# of reads supporting the variant allele."',
            )
            filtered_vcf.write(new_DV)
        elif line.startswith('##') and 'ID=CN' in line:
            continue
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if int(columns[headers.index('QUAL')]) >= qual and (columns[headers.index('FILTER')] == 'PASS'):
                if 'DUP' in columns[headers.index('ALT')]:
                    Format = columns[headers.index('FORMAT')].replace('DP', 'DR').replace('AD', 'DV')
                    Format_info = re.split(':|,', columns[headers.index(columnid)])
                    del Format_info[2]
                    filtered_vcf.write(
                        '{}\t{}\t{}\n'.format('\t'.join(columns[0:8]), Format, ':'.join(Format_info))
                    )
                else: 
                    Format = columns[headers.index('FORMAT')].replace('DP', 'DR').replace('AD', 'DV')
                    Format_info = re.split(':|,', columns[headers.index(columnid)])
                    del Format_info[1]
                    filtered_vcf.write(
                        '{}\t{}\t{}\n'.format('\t'.join(columns[0:8]), Format, ':'.join(Format_info))
                    )
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()


def reformat_sniffles(inp, out, sampleid, columnid):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('##') and 'ID=SEQ' in line:
            new_SEQ = line.replace(
                '##INFO=<ID=SEQ,Number=1,Type=String,Description="Extracted sequence from the best representative read.">',
                '##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Sequence of insertion">',
            )
            filtered_vcf.write(new_SEQ)
        elif line.startswith('##FILTER'):
            filtered_vcf.write('##FILTER=<ID=STRANDBIAS,Description="Strand is biased.">\n')
        elif line.startswith('##ALT') and 'TRA' in line:
            new_BND = line.replace('TRA,Description="Translocation"', 'BND,Description="Breakend"')
            filtered_vcf.write(new_BND)
        # elif line.startswith('##ALT') and 'INVDUP' in line:
        #     continue
        elif line.startswith('#CHROM'):
            headers = line.strip().replace(sampleid, columnid).split('\t')
            filtered_vcf.write(line.replace(sampleid, columnid))
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if columns[headers.index('FILTER')] == 'PASS' and 'IMPRECISE' not in line:
                if 'DEL' in columns[headers.index('INFO')]:
                    columns[headers.index('REF')] = 'N'
                    columns[headers.index('ALT')] = '<DEL>'
                    filtered_vcf.write('\t'.join(columns) + '\n')
                elif 'INS' in columns[headers.index('INFO')]:
                    # columns[headers.index('INFO')] += ';SVINSSEQ={}'.format(
                    #     columns[headers.index('ALT')]
                    # )
                    columns[headers.index('ALT')] += '.'
                    filtered_vcf.write('\t'.join(columns) + '\n')
                # elif 'INVDUP' in columns[headers.index('ALT')]:
                #     columns[headers.index('ALT')] = '<INV>'
                #     columns[headers.index('INFO')] = columns[headers.index('INFO')].replace('INVDUP', 'INV')
                #     filtered_vcf.write('\t'.join(columns) + '\n')
                else:
                    filtered_vcf.write(line)
        else:
            filtered_vcf.write(line)


def reformat_cutesv(inp, out):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(
                '##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Sequence of insertion">\n'
            )
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if columns[headers.index('FILTER')] == 'PASS' and 'IMPRECISE' not in line:
                if 'DEL' in columns[headers.index('INFO')]:
                    columns[headers.index('REF')] = 'N'
                    columns[headers.index('ALT')] = '<DEL>'
                    filtered_vcf.write('\t'.join(columns) + '\n')
                elif 'INS' in columns[headers.index('INFO')]:
                    columns[headers.index('REF')] = 'N'
                    columns[headers.index('INFO')] += ';SVINSSEQ={}'.format(
                        columns[headers.index('ALT')]
                    )
                    columns[headers.index('ALT')] = '<INS>'
                    filtered_vcf.write('\t'.join(columns) + '\n')
                else:
                    columns[headers.index('REF')] = 'N'
                    filtered_vcf.write('\t'.join(columns) + '\n')
        else:
            filtered_vcf.write(line)
