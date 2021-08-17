#! /usr/bin/python

'''
This pipeline computes somatic variants from WGS long read tumor-normal paired
data.py

This pipeline trims long-reads with porechop, aligns with minimap2,
performs structural variant calling with Sniffles, CuteSV, SVIM and
nanomonSV, and single nucleotide variants with Clai3, NanoCaller and
PEPPER.
The variants are then combined into one big file and annotated with VEP.py

Multiple filtering options are available. To see them type --help.

@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
'''

import datetime
import logging
import multiprocessing as mp
import os
import shutil
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from scripts.__version__ import version
from scripts.common import *
from scripts.reformat import *


def main(FQ_NORMAL, FQ_TUMOR, SAMPLEID, GENOME_REF, THREADS, STEPS, ASSEMBLY):

    logging.basicConfig(
        format='[%(asctime)s] - [%(levelname)s] - %(message)s',
        datefmt='%d-%b-%y %H:%M:%S',
        level=logging.DEBUG,
        filename=SAMPLEID + '.log',
    )
    logger = logging.getLogger(SAMPLEID)

    start_pipeline_time = datetime.datetime.now()

    logger.info(
        '###############################################################################################'
    )
    logger.info(
        '############# Starting WGS Long Read Somatic Pipeline: {} #############'.format(
            start_pipeline_time
        )
    )
    logger.info(
        '###############################################################################################'
    )
    logger.info('Pipeline version: {}'.format(version))
    logger.info(
        'Processing Normal FASTQ {}; and Tumor FASTQ {} with Sample ID {} '
        'using reference genome {}.'.format(FQ_NORMAL, FQ_TUMOR, SAMPLEID, ASSEMBLY)
    )

    sample_normal = SAMPLEID + '_Normal'
    sample_tumor = SAMPLEID + '_Tumor'

    os.makedirs('workdir', exist_ok=True)
    os.chdir('workdir')

    if 'mapping' in STEPS:
        SAM_THREADS = max(int(THREADS / 2), 1)

        start_map_time = datetime.datetime.now()
        logger.info('Starting trimming and mapping steps: {}'.format(start_map_time))

        # TRIMMING To-Do
        logger.info('Started trimming')

        # MAPPING
        logger.info('Starting alignment.')

        cmd = '{} -ax map-ont -t {} --MD {}.mmi {} -o {}.sam'.format(
            MINIMAP, THREADS, GENOME_REF, FQ_TUMOR, sample_tumor
        )

        p1 = exec_command(cmd, detach=True)

        cmd = '{} -ax map-ont -t {} --MD {}.mmi {} -o {}.sam'.format(
            MINIMAP, THREADS, GENOME_REF, FQ_NORMAL, sample_normal
        )

        p2 = exec_command(cmd, detach=True)

        p1.wait()
        p2.wait()

        cmd = '{} sort -@ {} -o {}.bam {}.sam'.format(
            SAMTOOLS, SAM_THREADS, sample_tumor, sample_tumor
        )
        p1 = exec_command(cmd, detach=True)

        cmd = '{} sort -@ {} -o {}.bam {}.sam'.format(
            SAMTOOLS, SAM_THREADS, sample_normal, sample_normal
        )
        p2 = exec_command(cmd, detach=True)

        p1.wait()
        p2.wait()

        cmd = '{} index {}.bam'.format(SAMTOOLS, sample_tumor)
        p1 = exec_command(cmd, detach = True)

        cmd = '{} index {}.bam'.format(SAMTOOLS, sample_normal)  
        p2 = exec_command(cmd, detach = True)

        p1.wait()
        p2.wait()


        end_map_time = datetime.datetime.now()
        total_map_time = end_map_time - start_map_time
        logger.info('Total trimming and mapping execution time: {}'.format(total_map_time))

    if 'variant' in STEPS:

        start_variant_time = datetime.datetime.now()
        logger.info('Starting variant calling: {}'.format(start_variant_time))
        logger.info('Variant calling with nanomonSV')

        cmd = '{} parse {}.bam nanomon_vc/Tumor'.format(NANOMON, sample_tumor)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} parse {}.bam nanomon_vc/Normal'.format(NANOMON, sample_normal)
        p2 = exec_command(cmd, detach=True)

        p1.wait()
        p2.wait()

        cmd = '{} get nanomon_vc/Tumor {}.bam {} --use_racon --control_prefix nanomon_vc/Normal --control_bam {}.bam'.format(
            NANOMON, sample_tumor, GENOME_REF, sample_normal
        )
        p3 = exec_command(cmd, detach=True)

        logger.info('Variant calling with Sniffles')

        cmd = '{} -s 2 -t {} --genotype --report_BND -m {}.bam -v {}_sniffles.vcf'.format(
            SNIFFLES, THREADS, sample_tumor, sample_tumor
        )
        p4 = exec_command(cmd, detach=True)

        cmd = '{} -s 2 -t {} --genotype --report_BND -m {}.bam -v {}_sniffles.vcf'.format(
            SNIFFLES, THREADS, sample_normal, sample_normal
        )
        p5 = exec_command(cmd, detach=True)

        logger.info('Variant calling with SVIM')

        # Call variants with SVIM for Tumor sample

        cmd = '{} --sample SVIM_Tumor --tandem_duplications_as_insertions --interspersed_duplications_as_insertions svim_tumor/ {}.bam {}'.format(
            SVIM, sample_tumor, GENOME_REF
        )
        p6 = exec_command(cmd, detach=True)

        # Call variants with SVIM for Normal sample

        cmd = '{} --sample SVIM_Normal --tandem_duplications_as_insertions --interspersed_duplications_as_insertions svim_normal/ {}.bam {}'.format(
            SVIM, sample_normal, GENOME_REF
        )
        p7 = exec_command(cmd, detach=True)

        logger.info('Variant calling with CuteSV')

        # Call variants with cuteSV for Tumor sample

        cmd = (
            '{} -t {} -S CUTESV_Tumor -s 2 --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 '
            '{}.bam {} CUTESV_Tumor.vcf .'.format(CUTESV, THREADS, sample_tumor, GENOME_REF)
        )
        p8 = exec_command(cmd, detach=True)

        # Call variants with cuteSV for Normal sample

        cmd = (
            '{} -t {} -S CUTESV_Normal -s 2 --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 '
            '{}.bam {} CUTESV_Normal.vcf .'.format(CUTESV, THREADS, sample_normal, GENOME_REF)
        )
        p9 = exec_command(cmd, detach=True)

        p3.wait()
        p4.wait()
        p5.wait()
        p6.wait()
        p7.wait()
        p8.wait()
        p9.wait()

        end_variant_time = datetime.datetime.now()
        total_variant_time = end_variant_time - start_variant_time
        logger.info('Total variant calling time: {}'.format(total_variant_time))

    if 'filter' in STEPS:

        start_filter_time = datetime.datetime.now()
        logger.info('Starting filtering and reformatting: {}'.format(start_filter_time))

        # Reformat SVIM VCF to follow Sniffles format and filter on QUAL

        p1 = mp.Process(
            target=reformat_svim,
            args=('svim_tumor/variants.vcf', 'tmp_svim_tumor.vcf', 'SVIM_Tumor', 0),
        )
        p2 = mp.Process(
            target=reformat_svim,
            args=('svim_normal/variants.vcf', 'tmp_svim_normal.vcf', 'SVIM_Normal', 0),
        )

        p1.start()
        p2.start()

        p1.join()
        p2.join()

        # Reformat nanomonsv VCF to follow Sniffles format

        p3 = mp.Process(
            target=reformat_nanomonsv,
            args=('nanomon_vc/Tumor.nanomonsv.result.vcf', 'tmp_nanomonsv.vcf'),
        )
        p3.start()
        p3.join()

        # Filter CuteSV VCF

        # p4 = mp.Process(target=, args=)
        # p4.start()
        # p4.join()

        ## TO-DO:
        # Add per-caller merge and variant filtering.
        # Add ensemble merge and filter based on number of callers.

        end_filter_time = datetime.datetime.now()
        total_filter_time = end_filter_time - start_filter_time
        logger.info('Total filtering and reformatting time: {}'.format(total_filter_time))

    if 'annotation' in STEPS:

        start_annotation_time = datetime.datetime.now()
        logger.info('Starting annotation: {}'.format(start_annotation_time))

        try:
            if SNPEFFDB == 'hg38':

                cmd = '{} -Xmx16g GRCh38.86 -csvStats combined_calls_snpEff.csv -v combined_calls.vcf > annotated_{}.vcf'.format(
                    SNPEFF, SNPEFFDB
                )
                exec_command(cmd)

            elif SNPEFFDB == 'hg19':

                cmd = '{} -Xmx16g GRCh37.75 -csvStats combined_calls_snpEff.csv -v combined_calls.vcf > annotated_{}.vcf'.format(
                    SNPEFF, SNPEFFDB
                )
                exec_command(cmd)

            else:
                raise ValueError(
                    "Wrong database specified, please specify one of these two: hg38 or hg19."
                )

        except ValueError as ve:
            print(ve)
            logger.error(ve)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('FASTQ_NORMAL', help='FASTQ file for Normal sample.')
    parser.add_argument('FASTQ_TUMOR', help='FASTQ file for Tumor sample')
    parser.add_argument(
        '--genome', type=str, required=True, help='Path to the reference genome FASTA file.'
    )
    parser.add_argument(
        '--sample',
        type=str,
        help='Name of the sample/experiment. Default is sample',
        default='sample',
    )
    parser.add_argument(
        '--outdir',
        type=str,
        required=True,
        help='Path to the output folder where output files will be placed',
    )
    parser.add_argument(
        '--snpeff-db',
        type=str,
        default='hg38',
        required=False,
        help='Genome assembly version to be used in snpEff (default: hg38)',
    )
    parser.add_argument(
        '--threads',
        help='Number of threads to use in the parallel steps',
        type=int,
        default=10,
        required=False,
    )
    parser.add_argument(
        '--steps',
        nargs='+',
        default=['mapping', 'variant', 'filter', 'annotation'],
        help='Steps to perform in the pipeline',
        choices=['mapping', 'variant', 'filter', 'annotation'],
    )
    # parser.add_argument('--keep-intermediate', default=False, action='store_true', required=False,
    #                     help='Do not remove temporary files')

    # Parse arguments
    args = parser.parse_args()
    DIR = args.outdir
    FQ_NORMAL = os.path.abspath(args.FASTQ_NORMAL)
    FQ_TUMOR = os.path.abspath(args.FASTQ_TUMOR)
    SAMPLEID = args.sample
    GENOME_REF = os.path.abspath(args.genome)
    THREADS = int(args.threads)
    STEPS = args.steps
    ASSEMBLY = args.snpeff_db

    # Move to output dir
    os.makedirs(os.path.abspath(DIR), exist_ok=True)
    os.chdir(os.path.abspath(DIR))

    main(FQ_NORMAL, FQ_TUMOR, SAMPLEID, GENOME_REF, THREADS, STEPS, ASSEMBLY)
