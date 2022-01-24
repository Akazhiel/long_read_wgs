"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
"""

from scripts.tools import *
from scripts.common import exec_command


def merge_variants(inp, out, distance):
    cmd = '{} --vcf {} --overlap 0.5 --ins_distance {} --bnd_distance {} > {}'.format(
        SVDB, ' '.join(inp), distance, distance, out
    )
    exec_command(cmd)
