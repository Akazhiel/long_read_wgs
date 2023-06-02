"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
"""

from scripts.tools import *
from scripts.common import exec_command


# def merge_variants(inp, out, distance):
#     cmd = '{} --vcf {} --overlap 0.95 --ins_distance 1 --bnd_distance {} > {}'.format(
#         SVDB, ' '.join(inp), distance, out
#     )
#     exec_command(cmd)

def merge_variants(inp, out, distance):
    with open("merge_vcfs.txt", "w") as handle:
        for s in inp:
            handle.write(s + "\n")
    cmd = '{} merge_vcfs.txt {} 1 1 1 -1 -1 {}'.format(
        SURVIRVOR, distance, out
    )
    exec_command(cmd)

    cmd = 'sed -i "s/:DR:/:AD:/" {}'.format(out)

    exec_command(cmd)

    cmd = 'sed -i "s/ID=DR/ID=AD/" {}'.format(out)

    exec_command(cmd)

def add_ins_sequence(merged_vcf, source_vcf, output_vcf):
    # Sort and index the merged vcf
    cmd = '{} sort -O z -o {}.gz {} && {} -s 1 -b 2 -e 2 {}.gz'.format(BCFTOOLS, merged_vcf, merged_vcf, TABIX, merged_vcf)
    sort1 = exec_command(cmd, detach = True)
    
    # Sort and index the source vcf
    cmd = '{} sort -O z -o {}.gz {} && {} -p vcf {}.gz'.format(BCFTOOLS, source_vcf, source_vcf, TABIX, source_vcf)
    sort2 = exec_command(cmd, detach = True)

    sort1.wait()
    sort2.wait()

    # Add insertion sequences to the merged vcf
    cmd = '{} annotate -a {}.gz -c INFO/SVINSSEQ -o {} {}.gz'.format(BCFTOOLS, source_vcf, output_vcf, merged_vcf)
    exec_command(cmd)