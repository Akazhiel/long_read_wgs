"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com> 
"""

import re
import logging

def filter_nanomonsv(inp, out):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('##') and 'ID=TR' in line:
            filtered_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            new_DR = line.replace('ID=TR,Number=1,Type=Integer,Description="The number of reads around the breakpoints"',
                        'ID=DR,Number=1,Type=Integer,Description="# of reads supporting the reference allele."')
            filtered_vcf.write(new_DR)
        elif line.startswith('##') and 'ID=VR' in line:
            new_DV = line.replace('ID=VR,Number=1,Type=Integer,Description="The number of variant supporting reads determined in the validation realignment step"',
                        'ID=DV,Number=1,Type=Integer,Description="# of reads supporting the variant allele."')
            filtered_vcf.write(new_DV)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            Format = columns[headers.index('FORMAT')].replace('TR', 'DR').replace('VR', 'DV')
            Format = 'GT:' + Format
            Normal = './.:' + columns[headers.index('CONTROL')]
            Tumor = './.:' + columns[headers.index('TUMOR')]
            filtered_vcf.write('{}\t{}\t{}\t{}\n'.format('\t'.join(columns[0:8]), Format, Tumor, Normal))
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()

def reformat_svim(inp, out, sample):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('##') and 'ID=DP' in line:
            new_DR = line.replace('ID=DP,Number=1,Type=Integer,Description="Read depth"',
                        'ID=DR,Number=1,Type=Integer,Description="# reads supporting the reference allele."')
            filtered_vcf.write(new_DR)
        elif line.startswith('##') and 'ID=AD' in line:
            new_DV = line.replace('ID=AD,Number=R,Type=Integer,Description="Read depth for each allele"',
                        'ID=DV,Number=R,Type=Integer,Description="# of reads supporting the variant allele."')
            filtered_vcf.write(new_DV)
        elif line.startswith('##') and 'ID=CN' in line:
            continue
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            Format = columns[headers.index('FORMAT')].replace('DP', 'DR').replace('AD', 'DV')
            Tumor = re.split(':|,', columns[headers.index(sample)])
            del Tumor[1]
            filtered_vcf.write('{}\t{}\t{}\n'.format('\t'.join(columns[0:8]), Format, ':'.join(Tumor)))
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()