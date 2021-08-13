#! /usr/bin/python

"""
This pipeline computes somatic variants from WGS long read tumor-normal paired 
data.py

This pipeline trims long-reads with porechop, aligns with minimap2,
performs structural variant calling with Sniffles, CuteSV, SVIM and 
nanomonSV, and single nucleotide variants with Clai3, NanoCaller and 
PEPPER.
The variants are then combined into one big file and annotated with VEP.py

Multiple filtering options are available. To see them type --help.

@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com> 
"""

import datetime
import logging
import multiprocessing as mp
import os
import shutil

from scripts.__version__ import version
from scripts.common import *
from scripts.reformat import *

SAMPLEID = "TEST"


def main(SAMPLEID):

    logging.basicConfig(
        format="%(asctime)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.DEBUG,
        filename=SAMPLEID + ".log",
    )
    logger = logging.getLogger(SAMPLEID)

    start_pipeline_time = datetime.datetime.now()

    logger.info(
        "###############################################################################################"
    )
    logger.info(
        "############# Starting WGS Long Read Somatic Pipeline: {} #############".format(
            start_pipeline_time
        )
    )
    logger.info(
        "###############################################################################################"
    )
    logger.info("Pipeline version: {}".format(version))
    logger.info(
        "Processing Normal FASTQ {}; and Tumor FASTQ {} with Sample ID {} "
        "using reference genome {}.".format(
            SAMPLEID, SAMPLEID, SAMPLEID, SAMPLEID
        )
    )

    sample_normal = SAMPLEID + "_Normal"
    sample_tumor = SAMPLEID + "_Tumor"

    os.makedirs("workdir", exist_ok=True)

    if "mapping" in STEPS:
        SAM_THREADS = max(int(THREADS / 2), 1)

        start_map_time = datetime.datetime.now()
        logger.info(
            "Starting trimming and mapping steps: {}".format(start_map_time)
        )

        # TRIMMING To-Do
        logger.info("Started trimming")

        # MAPPING
        logger.info("Starting alignment.")

        cmd = "{} -ax map-ont -t {} --MD {} {} -o {}.sam".format(
            MINIMAP, THREADS, GENOME, FASTQ_TUMOR, sample_tumor
        )

        p1 = exec_command(cmd)

        cmd = "{} -ax map-ont -t {} --MD {} {} -o {}.sam".format(
            MINIMAP, THREADS, GENOME, FASTQ_NORMAL, sample_normal
        )

        p2 = exec_command(cmd)

        p1.wait()
        p2.wait()

        cmd = "{} sort -@ {} -o {}.bam {}.sam".format(
            SAMTOOLS, SAM_THREADS, sample_tumor, sample_tumor
        )
        p1 = exec_command(cmd)

        cmd = "{} sort -@ {} -o {}.bam {}.sam".format(
            SAMTOOLS, SAM_THREADS, sample_normal, sample_normal
        )
        p2 = exec_command(cmd)

        p1.wait()
        p2.wait()

        end_map_time = datetime.datetime.now()
        total_map_time = end_map_time - start_map_time
        logger.info(
            "Total trimming and mapping execution time: {}".format(
                total_map_time
            )
        )

    if "variant" in STEPS:

        start_variant_time = datetime.datetime.now()
        logger.info("Starting variant calling: {}".format(start_variant_time))
        logger.info("Variant calling with nanomonSV")

        cmd = "{} parse {}.bam output/tumor".format(NANOMON, sample_tumor)
        p1 = exec_command(cmd)

        cmd = "{} parse {}.bam output/normal".format(NANOMON, sample_normal)
        p2 = exec_command(cmd)

        p1.wait()
        p2.wait()

        cmd = "{} get output/tumor {}.bam {} --control_prefix output/normal --control_bam {}.bam".format(
            NANOMON, sample_tumor, GENOME, sample_normal
        )
        p3 = exec_command(cmd)

        logger.info("Variant calling with Sniffles")

        cmd = "{} -s 2 -t {} --genotype --report_BND -m {}.bam -v {}_sniffles.vcf".format(
            SNIFFLES, THREADS, sample_tumor, sample_tumor
        )
        p4 = exec_command(cmd)

        cmd = "{} -s 2 -t {} --genotype --report_BND -m {}.bam -v {}_sniffles.vcf".format(
            SNIFFLES, THREADS, sample_normal, sample_normal
        )
        p5 = exec_command(cmd)

        logger.info("Variant calling with SVIM")

        cmd = "{} -t {} --sample SVIM_Tumor --tandem_duplications_as_insertions --interspersed_duplications_as_insertions {} {}.bam {}".format(
            SVIM, THREADS, ".", sample_tumor, GENOME
        )
        p6 = exec_command(cmd)

        cmd = "{} -t {} --sample SVIM_Normal --tandem_duplications_as_insertions --interspersed_duplications_as_insertions {} {}.bam {}".format(
            SVIM, THREADS, ".", sample_normal, GENOME
        )
        p7 = exec_command(cmd)

        logger.info("Variant calling with CuteSV")

        cmd = (
            "{} -t {} -S CUTESV_Tumor -s 2 --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 "
            "{} {} {} {}".format(CUTESV, THREADS)
        )


if __name__ == "__main__":
    main(SAMPLEID)
