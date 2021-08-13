"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com> 
"""

import subprocess
import sys
from scripts.tools import *
import os
import logging


def exec_command(cmd, detach=False):
    logger = logging.getLogger()
    logger.info(cmd)
    if detach:
        return subprocess.Popen(cmd, shell=True)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = p.communicate()
    if p.returncode != 0:
        for line in output.decode("utf-8").split("\n") if output else "":
            logger.error(line.rstrip())
        for line in error.decode("utf-8").split("\n") if error else "":
            logger.error(line.rstrip())
        sys.exit(-1)
