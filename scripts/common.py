'''
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
'''

import logging
import subprocess
import sys
import pandas as pd

from scripts.tools import *

class Variant:
    def __init__(self):
        self.chrom      = None
        self.start      = None
        self.end        = None
        self.alt        = None
        self.epitope    = None
        self.gene       = None
        self.transcript = None
        self.hgvsg      = None
        self.hgvsc      = None
        self.priority   = None
        self.mut        = None
        self.wt         = None
        self.frameshift = None
    
    @property
    def key(self):
        return f'{self.chrom}:{self.start}_{self.end}_{self.alt}'

def exec_command(cmd, detach=False):
    logger = logging.getLogger()
    logger.info(cmd)
    if detach:
        return subprocess.Popen(cmd, shell=True)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = p.communicate()
    if p.returncode != 0:
        for line in output.decode('utf-8').split('\n') if output else '':
            logger.error(line.rstrip())
        for line in error.decode('utf-8').split('\n') if error else '':
            logger.error(line.rstrip())
        sys.exit(-1)

def filter_annotsv(annot_sv_df):
    variants_df = pd.read_csv(annot_sv_df, sep = '\t').iloc[:, :35]
    variants_df[['Location1', 'Location2']] = variants_df['Location'].str.split('-', expand=True)
    variants_df = variants_df[variants_df.Location != 'txStart-txEnd']
    variants_df = variants_df[~variants_df.apply(lambda x: x['Location1'] == x['Location2'] and 'intron' in x['Location'], axis = 1)]
    variants_df = variants_df[~variants_df.apply(lambda x: x['SV_type'] == 'DEL' and 'txStart' in x['Location1'] , axis = 1)]
    variants_df['Priority'] = variants_df.apply(lambda x: 2 if (x['SV_type'] == 'DEL' and 'txEnd' in x['Location2']) else 3 if (x['SV_type'] not in ['INS', 'DEL']) else 1, axis = 1)
    return variants_df