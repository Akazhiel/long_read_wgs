'''
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
'''

import logging
import subprocess
import sys

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
        self.wt         = None
        self.mut        = None
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
