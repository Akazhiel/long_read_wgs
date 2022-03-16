"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
"""

from scripts.tools import *
from scripts.common import exec_command


def merge_variants(inp, out, distance):
    cmd = '{} --vcf {} --overlap 0 --ins_distance 5 --bnd_distance {} > {}'.format(
        SVDB, ' '.join(inp), distance, out
    )
    exec_command(cmd)
