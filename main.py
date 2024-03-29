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
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from pyensembl import EnsemblRelease

from scripts.__version__ import version
from scripts.common import exec_command
# from scripts.epitope import create_epitope
from scripts.filters import filter_callers, filter_somatic, prioritize_variants
# from scripts.hgvs_notations import add_variant_hgvs
from scripts.reformat import *
from scripts.tools import *
from scripts.vcfmerge import add_ins_sequence, merge_variants


def main(FQ_NORMAL, FQ_TUMOR, SAMPLEID, GENOME_REF, THREADS, STEPS, NUM_CALLERS, WINDOW, ENSEMBL_VERSION):

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
        'using reference genome GRCh38.'.format(FQ_NORMAL, FQ_TUMOR, SAMPLEID)
    )

    sample_normal = SAMPLEID + '_Normal'
    sample_tumor = SAMPLEID + '_Tumor'

    os.makedirs('workdir', exist_ok=True)
    os.chdir('workdir')

    if 'mapping' in STEPS:

        start_map_time = datetime.datetime.now()
        logger.info(
            'Starting trimming and mapping steps: {}'.format(start_map_time))

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
        logger.info(
            'Total trimming and mapping execution time: {}'.format(total_map_time))

    if 'variant' in STEPS:

        start_variant_time = datetime.datetime.now()
        logger.info('Starting variant calling: {}'.format(start_variant_time))

        logger.info('Variant calling with SVIM')

        # Call variants with SVIM for Tumor sample

        cmd = '{} --sample SVIM_Tumor --min_sv_size 30 --tandem_duplications_as_insertions --interspersed_duplications_as_insertions --partition_max_distance 5 --cluster_max_distance 0.8 --max_consensus_length 10000000 svim_tumor/ {}.bam {}'.format(
            SVIM, sample_tumor, GENOME_REF
        )
        p4 = exec_command(cmd, detach=True)

        # Call variants with SVIM for Normal sample

        cmd = '{} --sample SVIM_Normal --min_sv_size 10 --tandem_duplications_as_insertions --interspersed_duplications_as_insertions --partition_max_distance 5 --cluster_max_distance 0.8 --max_consensus_length 10000000 svim_normal/ {}.bam {}'.format(
            SVIM, sample_normal, GENOME_REF
        )
        p5 = exec_command(cmd, detach=True)

        logger.info('Variant calling with CuteSV')

        # Call variants with cuteSV for Tumor sample

        os.makedirs('cutesv_tumor')

        cmd = (
            '{} -t {} -S CUTESV_Tumor -l 30 -s 7 -L -1 -md 5 -mi 5 --genotype --max_cluster_bias_INS 5 --diff_ratio_merging_INS 0.2 --max_cluster_bias_DEL 5 '
            '--diff_ratio_merging_DEL 0.2 --remain_reads_ratio 0.8 {}.bam {} CUTESV_Tumor.vcf cutesv_tumor/'.format(
                CUTESV, THREADS, sample_tumor, GENOME_REF
            )
        )
        p6 = exec_command(cmd, detach=True)

        # Call variants with cuteSV for Normal sample

        os.makedirs('cutesv_normal')

        cmd = (
            '{} -t {} -S CUTESV_Normal -l 10 -s 2 -L -1 -md 5 -mi 5 --genotype --max_cluster_bias_INS 5 --diff_ratio_merging_INS 0.2 --max_cluster_bias_DEL 5 '
            '--diff_ratio_merging_DEL 0.2 --remain_reads_ratio 0.8 {}.bam {} CUTESV_Normal.vcf cutesv_normal/'.format(
                CUTESV, THREADS, sample_normal, GENOME_REF
            )
        )
        p7 = exec_command(cmd, detach=True)

        # Call variants with SNIFFLES

        logger.info('Variant calling with Sniffles')

        cmd = '{} --sample-id SNIFFLES_Tumor --minsupport 7 --allow-overwrite --cluster-binsize 5 --cluster-merge-len 0 --quiet --reference {} -t {} --minsvlen 30 --mapq 20 --input {}.bam --vcf {}_sniffles.vcf'.format(
            SNIFFLES, GENOME_REF, THREADS, sample_tumor, sample_tumor
        )
        p8 = exec_command(cmd, detach=True)

        cmd = '{} --sample-id SNIFFLES_Normal --minsupport 2 --allow-overwrite --cluster-binsize 5 --cluster-merge-len 0 --no-consensus --no-qc --quiet --reference {} -t {} --minsvlen 10 --mapq 20 --input {}.bam --vcf {}_sniffles.vcf'.format(
            SNIFFLES, GENOME_REF, THREADS, sample_normal, sample_normal
        )
        p9 = exec_command(cmd, detach=True)
        
        # Wait for all variant calling processes to finish
        
        p4.wait()
        p5.wait()
        p6.wait()
        p7.wait()
        p8.wait()
        p9.wait()

        end_variant_time = datetime.datetime.now()
        total_variant_time = end_variant_time - start_variant_time
        logger.info('Total variant calling time: {}'.format(
            total_variant_time))

    if 'filter' in STEPS:

        start_filter_time = datetime.datetime.now()
        logger.info('Starting filtering and reformatting: {}'.format(
            start_filter_time))

        # Reformat SVIM VCF to follow Sniffles format and filter on QUAL

        re_svimt = mp.Process(
            target=reformat_svim,
            args=('svim_tumor/variants.vcf',
                  'tmp_svim_tumor.vcf', 'SVIM_Tumor', 10),
        )
        re_svimn = mp.Process(
            target=reformat_svim,
            args=('svim_normal/variants.vcf',
                  'tmp_svim_normal.vcf', 'SVIM_Normal', 1),
        )

        re_svimt.start()
        re_svimn.start()

        re_svimt.join()
        re_svimn.join()

        # Filter variants in tumor vcf to keep only precise variants.

        cmd = f'{BCFTOOLS} view -i "INFO/STD_POS=\'.\'" tmp_svim_tumor.vcf > precise_svim_tumor.vcf'
        p1 = exec_command(cmd, detach=True)

        p1.wait()

        # Merge individual SVIM calls and filter for somatic variants

        merge_svim = mp.Process(
            target=merge_variants,
            args=(['precise_svim_tumor.vcf', 'tmp_svim_normal.vcf'],
                  'svim_combined_calls.vcf', 50)
        )

        merge_svim.start()

        # Reformat Sniffles to add INS Seq on INFO field

        cmd = f'{BCFTOOLS} view -i "STDEV_POS=0.0" {sample_tumor}_sniffles.vcf > precise_sniffles_tumor.vcf'
        p2 = exec_command(cmd, detach=True)

        # p5.wait()
        p2.wait()

        re_snifflesT = mp.Process(
            target=reformat_sniffles,
            args=(
                '{}_sniffles.vcf'.format(sample_normal),
                'tmp_sniffles_normal.vcf',
                'SNIFFLES_Normal'
            ),
        )
        re_snifflesN = mp.Process(
            target=reformat_sniffles,
            args=(
                'precise_sniffles_tumor.vcf',
                'tmp_sniffles_tumor.vcf',
                'SNIFFLES_Tumor'
            ),
        )

        re_snifflesT.start()
        re_snifflesN.start()

        re_snifflesT.join()
        re_snifflesN.join()

        # Merge individual SNIFFLES calls and filter for somatic variants

        merge_sniffles = mp.Process(
            target=merge_variants,
            args=(['tmp_sniffles_tumor.vcf', 'tmp_sniffles_normal.vcf'],
                  'sniffles_combined_calls.vcf', 50)
        )

        merge_sniffles.start()

        # Reformat CuteSV

        cmd = f'{BCFTOOLS} view -i "CIPOS=0" CUTESV_Tumor.vcf > precise_cutesv_tumor.vcf'
        p3 = exec_command(cmd, detach=True)

        # p9.wait()
        p3.wait()

        re_cutesvT = mp.Process(
            target=reformat_cutesv,
            args=('CUTESV_Normal.vcf', 
                  'tmp_cutesv_normal.vcf',
                  'CUTESV_Normal'
                  )
        )

        re_cutesvN = mp.Process(
            target=reformat_cutesv,
            args=('precise_cutesv_tumor.vcf', 
                  'tmp_cutesv_tumor.vcf',
                  'CUTESV_Tumor'
                  )
        )

        re_cutesvT.start()
        re_cutesvN.start()

        re_cutesvT.join()
        re_cutesvN.join()

        # Merge individual CUTESV calls and filter for somatic variants

        merge_cutesv = mp.Process(
            target=merge_variants,
            args=(['tmp_cutesv_tumor.vcf', 'tmp_cutesv_normal.vcf'],
                  'cutesv_combined_calls.vcf', 50)
        )

        merge_cutesv.start()

        # Wait for merges to finish

        merge_svim.join()
        merge_sniffles.join()
        merge_cutesv.join()

        # Add insertion sequences to merged VCFs

        add_ins_sequence('svim_combined_calls.vcf', 'precise_svim_tumor.vcf', 'svim_combined.vcf')
        add_ins_sequence('sniffles_combined_calls.vcf', 'tmp_sniffles_tumor.vcf', 'sniffles_combined.vcf')
        add_ins_sequence('cutesv_combined_calls.vcf', 'tmp_cutesv_tumor.vcf', 'cutesv_combined.vcf')

        # Filter somatic calls for each caller

        filter_somatic('svim_combined.vcf',
                       'svim_combined_filtered.vcf',
                       'SVIM'
                       )

        filter_somatic('sniffles_combined.vcf',
                       'sniffles_combined_filtered.vcf',
                       'SNIFFLES'
                       )

        filter_somatic('cutesv_combined.vcf',
                       'cutesv_combined_filtered.vcf',
                       'CUTESV'
                       )

        # Merge the combined calls from the different callers and filter based on number of callers

        merge_variants(
            [
                'sniffles_combined_filtered.vcf',
                'svim_combined_filtered.vcf',
                'cutesv_combined_filtered.vcf',
            ],
            'combined_calls.vcf',
            WINDOW,
        )

        filter_callers('combined_calls.vcf',
                       'combined_calls_filtered.vcf', NUM_CALLERS)

        end_filter_time = datetime.datetime.now()
        total_filter_time = end_filter_time - start_filter_time
        logger.info('Total filtering and reformatting time: {}'.format(
            total_filter_time))

    if 'annotation' in STEPS:

        start_annotation_time = datetime.datetime.now()
        logger.info('Starting annotation: {}'.format(start_annotation_time))

        # Annotate variants using AnnotSV with latest ENSEMBL release

        cmd = f'{ANNOTSV} -tx ENSEMBL -SVinputFile combined_calls_filtered.vcf -SVminSize 30 -outputFile annotsv_ensembl.tsv -outputDir . -annotationMode split -genomeBuild GRCh38'
        exec_command(cmd)

        end_annotation_time = datetime.datetime.now()
        total_annotation_time = end_annotation_time - start_annotation_time
        logger.info('Total annotation time: {}'.format(total_annotation_time))

        # Prirotize variants according to breakpoints.
    if 'epitope' in STEPS:

        ensembl_data = EnsemblRelease(ENSEMBL_VERSION)

        annotsv_prio = prioritize_variants('annotsv_ensembl.tsv')

        variants = add_variant_hgvs(annotsv_prio)

        for variant in variants:
            create_epitope(variant)

    end_pipeline_time = datetime.datetime.now()
    total_pipeline_time = end_pipeline_time - start_pipeline_time
    logger.info('Total pipeline execution time: {}'.format(
        total_pipeline_time))

    logger.info('COMPLETED!')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
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
        default=['mapping', 'variant', 'filter', 'annotation', 'epitope'],
        help='Steps to perform in the pipeline. List of choices: {%(choices)s}',
        choices=['mapping', 'variant', 'filter', 'annotation', 'epitope'],
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
    parser.add_argument(
        '--ensembl-version',
        metavar='\b',
        help='Ensembl version number that was used to annotate the variants with AnnotSV. (Default: %(default)s)',
        type=int,
        default=107,
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
    NUM_CALLERS = args.num_callers
    WINDOW = args.window
    ENSEMBL_VERSION = args.ensembl_version

    # Move to output dir
    os.makedirs(os.path.abspath(DIR), exist_ok=True)
    os.chdir(os.path.abspath(DIR))

    main(FQ_NORMAL, FQ_TUMOR, SAMPLEID, GENOME_REF, THREADS,
         STEPS, NUM_CALLERS, WINDOW, ENSEMBL_VERSION)
