#! /usr/bin/python

'''
This pipeline computes somatic variants from WGS long read tumor-normal paired
data.py

This pipeline trims long-reads with porechop, aligns with minimap2,
performs structural variant calling with Sniffles, CuteSV, SVIM and
nanomonSV, and single nucleotide variants with Clair3, NanoCaller and
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
import time
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from scripts.__version__ import version
from scripts.common import *
from scripts.filters import filter_somatic, filter_callers
from scripts.reformat import *
from scripts.vcfmerge import merge_variants


def main(FQ_NORMAL, FQ_TUMOR, SAMPLEID, GENOME_REF, THREADS, STEPS, SNPEFFDB, NUM_CALLERS, WINDOW):

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
        'using reference genome {}.'.format(FQ_NORMAL, FQ_TUMOR, SAMPLEID, SNPEFFDB)
    )

    sample_normal = SAMPLEID + '_Normal'
    sample_tumor = SAMPLEID + '_Tumor'

    os.makedirs('workdir', exist_ok=True)
    os.chdir('workdir')

    if 'mapping' in STEPS:

        start_map_time = datetime.datetime.now()
        logger.info('Starting trimming and mapping steps: {}'.format(start_map_time))

        # TRIMMING To-Do
        logger.info('Started trimming')

        # MAPPING
        logger.info('Starting alignment.')

        SAMTOOLS_THREADS = max(int(THREADS / 2), 1)

        cmd = '{} -ax map-ont -t {} -K4G --MD {}.mmi {} | {} sort -m 2G --threads {} -o {}.bam'.format(
            MINIMAP, THREADS, GENOME_REF, FQ_TUMOR, SAMTOOLS, SAMTOOLS_THREADS, sample_tumor
        )

        p1 = exec_command(cmd, detach=True)

        cmd = '{} -ax map-ont -t {} --MD -K4G {}.mmi {} | {} sort -m 2G --threads {} -o {}.bam'.format(
            MINIMAP, THREADS, GENOME_REF, FQ_NORMAL, SAMTOOLS, SAMTOOLS_THREADS, sample_normal
        )

        p2 = exec_command(cmd, detach=True)

        p1.wait()
        p2.wait()

        cmd = '{} index {}.bam'.format(SAMTOOLS, sample_tumor)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} index {}.bam'.format(SAMTOOLS, sample_normal)
        p2 = exec_command(cmd, detach=True)

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

        time.sleep(1)  # Required so the parse step by nanomon is executed in parallel but doesn't break due to the same folder name

        cmd = '{} parse {}.bam nanomon_vc/Normal'.format(NANOMON, sample_normal)
        p2 = exec_command(cmd, detach=True)

        p1.wait()
        p2.wait()

        cmd = '{} get nanomon_vc/Tumor {}.bam {} --var_read_min_mapq 20 --use_racon --control_prefix nanomon_vc/Normal --control_bam {}.bam'.format(
            NANOMON, sample_tumor, GENOME_REF, sample_normal
        )
        p3 = exec_command(cmd, detach=True)

        logger.info('Variant calling with SVIM')

        # Call variants with SVIM for Tumor sample

        cmd = '{} --sample SVIM_Tumor --min_sv_size 30 --max_consensus_length 1000000 --insertion_sequences --tandem_duplications_as_insertions --interspersed_duplications_as_insertions --symbolic_alleles svim_tumor/ {}.bam {}'.format(
            SVIM, sample_tumor, GENOME_REF
        )
        p4 = exec_command(cmd, detach=True)

        # Call variants with SVIM for Normal sample

        cmd = '{} --sample SVIM_Normal --min_sv_size 30 --max_consensus_length 1000000 --insertion_sequences --tandem_duplications_as_insertions --interspersed_duplications_as_insertions --symbolic_alleles svim_normal/ {}.bam {}'.format(
            SVIM, sample_normal, GENOME_REF
        )
        p5 = exec_command(cmd, detach=True)

        logger.info('Variant calling with CuteSV')

        # Call variants with cuteSV for Tumor sample

        os.makedirs('cutesv_tumor')

        cmd = (
            '{} -t {} -S CUTESV_Tumor -s 1 -L -1 --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 '
            '--diff_ratio_merging_DEL 0.5 {}.bam {} CUTESV_Tumor.vcf cutesv_tumor/'.format(
                CUTESV, THREADS, sample_tumor, GENOME_REF
            )
        )
        p6 = exec_command(cmd, detach=True)

        # Call variants with cuteSV for Normal sample

        os.makedirs('cutesv_normal')

        cmd = (
            '{} -t {} -S CUTESV_Normal -s 1 -L -1 --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 '
            '--diff_ratio_merging_DEL 0.5 {}.bam {} CUTESV_Normal.vcf cutesv_normal/'.format(
                CUTESV, THREADS, sample_normal, GENOME_REF
            )
        )
        p7 = exec_command(cmd, detach=True)

        p3.wait()
        p4.wait()
        p5.wait()
        p6.wait()
        p7.wait()

        # Call variants with SNIFFLES

        logger.info('Variant calling with Sniffles')

        cmd = '{} --minsupport 1 --symbolic --reference {} -t {} --genotype --minsvlen 30 --mapq 20 --qc-stdev-abs-max 0 --cluster-merge-pos 50 --input {}.bam --vcf {}_sniffles.vcf'.format(
            SNIFFLES, GENOME_REF, THREADS, sample_tumor, sample_tumor
        )
        p8 = exec_command(cmd, detach=True)
        p8.wait()

        cmd = '{} --minsupport 1 --symbolic --reference {} -t {} --genotype --minsvlen 30 --mapq 20 --qc-stdev-abs-max 0 --cluster-merge-pos 50 --input {}.bam --vcf {}_sniffles.vcf'.format(
            SNIFFLES, GENOME_REF, THREADS, sample_normal, sample_normal
        )
        p9 = exec_command(cmd, detach=True)
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
            args=('svim_tumor/variants.vcf', 'tmp_svim_tumor.vcf', 'SVIM_Tumor', 10),
        )
        p2 = mp.Process(
            target=reformat_svim,
            args=('svim_normal/variants.vcf', 'tmp_svim_normal.vcf', 'SVIM_Normal', 10),
        )

        p1.start()
        p2.start()

        p1.join()
        p2.join()

        # Merge individual SVIM calls and filter for somatic variants

        merge_variants(
            ['tmp_svim_tumor.vcf', 'tmp_svim_normal.vcf'], 'svim_combined_calls.vcf', 100
        )
        filter_somatic('svim_combined_calls.vcf', 'svim_combined_calls_filtered.vcf', 'SVIM')

        # Reformat nanomonsv VCF to follow Sniffles format

        genomesv2vcf_convert('nanomon_vc/Tumor.nanomonsv.result.txt', 'Tumor.nanomonsv.result.vcf', GENOME_REF)

        reformat_nanomonsv('Tumor.nanomonsv.result.vcf', 'tmp_nanomonsv.vcf')

        filter_somatic('tmp_nanomonsv.vcf', 'nanomonsv_filtered.vcf', 'NANOMON')

        # Reformat Sniffles to add INS Seq on INFO field

        p4 = mp.Process(
            target=reformat_sniffles,
            args=(
                '{}_sniffles.vcf'.format(sample_normal),
                'tmp_sniffles_normal.vcf',
                '{}.bam'.format(sample_normal),
                'SNIFFLES_Normal',
            ),
        )
        p5 = mp.Process(
            target=reformat_sniffles,
            args=(
                '{}_sniffles.vcf'.format(sample_tumor),
                'tmp_sniffles_tumor.vcf',
                '{}.bam'.format(sample_tumor),
                'SNIFFLES_Tumor',
            ),
        )

        p4.start()
        p5.start()

        p4.join()
        p5.join()

        # Merge individual SNIFFLES calls and filter for somatic variants

        merge_variants(
            ['tmp_sniffles_tumor.vcf', 'tmp_sniffles_normal.vcf'],
            'sniffles_combined_calls.vcf',
            100,
        )

        filter_somatic(
            'sniffles_combined_calls.vcf', 'sniffles_combined_calls_filtered.vcf', 'SNIFFLES'
        )

        # Reformat CuteSV

        p6 = mp.Process(target=reformat_cutesv, args=('CUTESV_Normal.vcf', 'tmp_cutesv_normal.vcf'))
        p7 = mp.Process(target=reformat_cutesv, args=('CUTESV_Tumor.vcf', 'tmp_cutesv_tumor.vcf'))

        p6.start()
        p7.start()

        p6.join()
        p7.join()

        # Merge individual CUTESV calls and filter for somatic variants

        merge_variants(
            ['tmp_cutesv_tumor.vcf', 'tmp_cutesv_normal.vcf'], 'cutesv_combined_calls.vcf', 100
        )

        filter_somatic('cutesv_combined_calls.vcf', 'cutesv_combined_calls_filtered.vcf', 'CUTESV')

        # Merge the combined calls from the different callers and filter based on number of callers

        merge_variants(
            [
                'nanomonsv_filtered.vcf',
                'svim_combined_calls_filtered.vcf',
                'sniffles_combined_calls_filtered.vcf',
                'cutesv_combined_calls_filtered.vcf',
            ],
            'combined_calls.vcf',
            WINDOW,
        )

        filter_callers('combined_calls.vcf', 'combined_calls_filtered.vcf', NUM_CALLERS)

        end_filter_time = datetime.datetime.now()
        total_filter_time = end_filter_time - start_filter_time
        logger.info('Total filtering and reformatting time: {}'.format(total_filter_time))

    if 'annotation' in STEPS:

        start_annotation_time = datetime.datetime.now()
        logger.info('Starting annotation: {}'.format(start_annotation_time))

        try:
            assert SNPEFFDB in [
                'hg38',
                'hg19',
            ], 'Database error: Wrong database provided, got {}; expected one of these two: hg38 or hg19.'.format(
                SNPEFFDB
            )

            if SNPEFFDB == 'hg38':

                cmd = '{} -Xmx16g -csvStats combined_calls_snpEff.csv -v GRCh38.99 combined_calls_filtered.vcf > annotated_{}.vcf'.format(
                    SNPEFF, SNPEFFDB
                )
                exec_command(cmd)

            elif SNPEFFDB == 'hg19':

                cmd = '{} -Xmx16g -csvStats combined_calls_snpEff.csv -v GRCh37.75 combined_calls_filtered.vcf > annotated_{}.vcf'.format(
                    SNPEFF, SNPEFFDB
                )
                exec_command(cmd)

        except AssertionError as ae:
            print(ae)
            logger.error(ae)
            sys.exit(-1)

        end_annotation_time = datetime.datetime.now()
        total_annotation_time = end_annotation_time - start_annotation_time
        logger.info('Total annotation time: {}'.format(total_annotation_time))

    end_pipeline_time = datetime.datetime.now()
    total_pipeline_time = end_pipeline_time - start_pipeline_time
    logger.info('Total pipeline execution time: {}'.format(total_pipeline_time))

    logger.info('COMPLETED!')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('FASTQ_NORMAL', help='FASTQ file for Normal sample.')
    parser.add_argument('FASTQ_TUMOR', help='FASTQ file for Tumor sample')
    parser.add_argument(
        '--genome',
        metavar='\b',
        type=str,
        required=True,
        help='Path to the reference genome FASTA file.',
    )
    parser.add_argument(
        '--sample',
        metavar='\b',
        type=str,
        help='Name of the sample/experiment. Default is %(default)s',
        default='sample',
    )
    parser.add_argument(
        '-o',
        '--outdir',
        metavar='\b',
        type=str,
        required=True,
        help='Path to the output folder where output files will be placed',
    )
    parser.add_argument(
        '-db',
        '--snpeff-db',
        metavar='\b',
        type=str,
        default='hg38',
        required=False,
        help='Genome assembly version to be used in snpEff. (Default: %(default)s)',
    )
    parser.add_argument(
        '-t',
        '--threads',
        metavar='\b',
        help='Number of threads to use in the parallel steps',
        type=int,
        default=10,
        required=False,
    )
    parser.add_argument(
        '--steps',
        nargs='+',
        default=['mapping', 'variant', 'filter', 'annotation'],
        help='Steps to perform in the pipeline. List of choices: {%(choices)s}',
        choices=['mapping', 'variant', 'filter', 'annotation'],
    )
    parser.add_argument(
        '-n',
        '--num_callers',
        metavar='\b',
        help='Filter for number of SV callers required. (Default: %(default)s)',
        type=int,
        default=2,
        required=False,
    )
    parser.add_argument(
        '-w',
        '--window',
        metavar='\b',
        help='Window threshold that is allowed to cluster two variants as the same. (Default: %(default)s)',
        type=int,
        default=50,
        required=False,
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
    SNPEFFDB = args.snpeff_db
    NUM_CALLERS = args.num_callers
    WINDOW = args.window

    # Move to output dir
    os.makedirs(os.path.abspath(DIR), exist_ok=True)
    os.chdir(os.path.abspath(DIR))

    main(FQ_NORMAL, FQ_TUMOR, SAMPLEID, GENOME_REF, THREADS, STEPS, SNPEFFDB, NUM_CALLERS, WINDOW)
