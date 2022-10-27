"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
"""

import re
import csv, datetime
import pysam


def reformat_nanomonsv(inp, out):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('##') and 'ID=TR' in line:
            filtered_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            new_DR = line.replace(
                'ID=TR,Number=1,Type=Integer,Description="The number of reads around the breakpoints"',
                'ID=DR,Number=1,Type=Integer,Description="# of reads supporting the reference allele."',
            )
            filtered_vcf.write(new_DR)
        elif line.startswith('##') and 'ID=VR' in line:
            new_DV = line.replace(
                'ID=VR,Number=1,Type=Integer,Description="The number of variant supporting reads determined in the validation realignment step"',
                'ID=DV,Number=1,Type=Integer,Description="# of reads supporting the variant allele."',
            )
            filtered_vcf.write(new_DV)
        elif line.startswith('#CHROM'):
            headers = (
                line.strip()
                .replace('TUMOR', 'NANOMON_Tumor')
                .replace('CONTROL', 'NANOMON_Normal')
                .split('\t')
            )
            filtered_vcf.write(
                line.replace('TUMOR', 'NANOMON_Tumor').replace('CONTROL', 'NANOMON_Normal')
            )
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if columns[headers.index('FILTER')] == 'PASS':
                columns[headers.index('REF')] = 'N'
                Format = columns[headers.index('FORMAT')].replace('TR', 'DR').replace('VR', 'DV')
                Format = 'GT:' + Format
                Normal_Format = columns[headers.index('NANOMON_Normal')].split(':')
                Normal_Format[0] = str(int(Normal_Format[0]) - int(Normal_Format[1]))
                Tumor_Format = columns[headers.index('NANOMON_Tumor')].split(':')
                Tumor_Format[0] = str(int(Tumor_Format[0]) - int(Tumor_Format[1]))
                Normal = './.:' + ':'.join(Normal_Format)
                Tumor = './.:' + ':'.join(Tumor_Format)
                filtered_vcf.write(
                    '{}\t{}\t{}\t{}\n'.format('\t'.join(columns[0:8]), Format, Tumor, Normal)
                )
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()


def reformat_svim(inp, out, columnid, qual):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('##') and 'ID=DP' in line:
            new_DR = line.replace(
                'ID=DP,Number=1,Type=Integer,Description="Read depth"',
                'ID=DR,Number=1,Type=Integer,Description="# reads supporting the reference allele."',
            )
            filtered_vcf.write(new_DR)
        elif line.startswith('##') and 'ID=AD' in line:
            new_DV = line.replace(
                'ID=AD,Number=R,Type=Integer,Description="Read depth for each allele"',
                'ID=DV,Number=1,Type=Integer,Description="# of reads supporting the variant allele."',
            )
            filtered_vcf.write(new_DV)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if int(columns[headers.index('QUAL')]) >= qual:
                Format = (
                        columns[headers.index('FORMAT')].replace('DP', 'DR').replace('AD', 'DV')
                    )
                Format_info = re.split(':|,', columns[headers.index(columnid)])
                del Format_info[1]
                if 'DUP:TANDEM' in columns[headers.index('ALT')]:
                    columns[headers.index('ALT')] = '<DUP>'
                    del Format[1]
                    del Format_info[1]
                    filtered_vcf.write(
                        '{}\t{}\t{}\n'.format(
                            '\t'.join(columns[0:8]), ':'.join(Format), ':'.join(Format_info)
                        )
                    )
                elif 'DUP:INT' in columns[headers.index('ALT')]:
                    columns[headers.index('ALT')] = '<DUP>'
                    filtered_vcf.write(
                        '{}\t{}\t{}\n'.format(
                            '\t'.join(columns[0:8]), Format, ':'.join(Format_info)
                        )
                    )
                elif 'DEL' in columns[headers.index('INFO')]:
                    columns[headers.index('REF')] = 'N'
                    columns[headers.index('ALT')] = '<DEL>'
                    columns[headers.index('POS')] = str(int(columns[headers.index('POS')]) + 1)
                    filtered_vcf.write(
                        '{}\t{}\t{}\n'.format(
                            '\t'.join(columns[0:8]), Format, ':'.join(Format_info)
                        )
                    )
                elif 'INS' in columns[headers.index('INFO')]:
                    columns[headers.index('INFO')] += ';SVINSSEQ={}'.format(columns[headers.index('ALT')])
                    columns[headers.index('ALT')] = '<INS>'
                    filtered_vcf.write(
                            '{}\t{}\t{}\n'.format(
                                '\t'.join(columns[0:8]), Format, ':'.join(Format_info)
                            )
                        )
                else:
                    filtered_vcf.write(
                            '{}\t{}\t{}\n'.format(
                                '\t'.join(columns[0:8]), Format, ':'.join(Format_info)
                            )
                        )

        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()


def reformat_sniffles(inp, out):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            new_SEQ = '##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Sequence of insertion">\n'
            filtered_vcf.write(new_SEQ)
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if 'DEL' in columns[headers.index('INFO')]:
                columns[headers.index('REF')] = 'N'
                columns[headers.index('ALT')] = '<DEL>'
                filtered_vcf.write('\t'.join(columns) + '\n')
            elif 'INS' in columns[headers.index('INFO')]:
                columns[headers.index('POS')] = str(int(columns[headers.index('POS')]) - 1)
                INFO = columns[headers.index('INFO')].split(';')
                pos_idx = [i for i, x in enumerate(INFO) if x.startswith('END')][0]
                INFO[pos_idx] = 'END=' + str(int(INFO[pos_idx].split('=')[1]) - 1)
                columns[headers.index('INFO')] = ';'.join(INFO)
                columns[headers.index('INFO')] += ';SVINSSEQ={}'.format(
                    columns[headers.index('ALT')]
                )
                columns[headers.index('ALT')] = '<INS>'
                filtered_vcf.write('\t'.join(columns) + '\n')
            else:
                filtered_vcf.write(line)
        else:
            filtered_vcf.write(line)


def reformat_cutesv(inp, out):
    vcf = open(inp, 'r')
    filtered_vcf = open(out, 'w')
    for line in vcf:
        if line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(
                '##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Sequence of insertion">\n'
            )
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if columns[headers.index('QUAL')] != '.':
                if 'DEL' in columns[headers.index('INFO')]:
                    columns[headers.index('REF')] = 'N'
                    columns[headers.index('ALT')] = '<DEL>'
                    columns[headers.index('POS')] = str(int(columns[headers.index('POS')]) + 1)
                    filtered_vcf.write('\t'.join(columns) + '\n')
                elif 'INS' in columns[headers.index('INFO')]:
                    columns[headers.index('POS')] = str(int(columns[headers.index('POS')]) - 1)
                    columns[headers.index('REF')] = 'N'
                    columns[headers.index('INFO')] += ';SVINSSEQ={}'.format(
                        columns[headers.index('ALT')]
                    )
                    columns[headers.index('ALT')] = '<INS>'
                    filtered_vcf.write('\t'.join(columns) + '\n')
                else:
                    columns[headers.index('REF')] = 'N'
                    filtered_vcf.write('\t'.join(columns) + '\n')
        else:
            filtered_vcf.write(line)

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))

def genomesv2vcf_convert(result_file, output_vcf, reference):

    today_str = datetime.datetime.today().strftime("%Y%m%d")

    header = '##fileformat=VCFv4.3\n'\
             f'##fileDate={today_str}\n'\
             f'##reference={reference}'

    ref_tb = pysam.FastaFile(reference)

    for (tchr, tlen) in zip(ref_tb.references, ref_tb.lengths):
        header = header + '\n' + f"##contig=<ID={tchr},length={tlen}>"

    header = header + '\n' + \
            '##FILTER=<ID=Duplicate_with_close_SV,Description="When multiple SVs that share breakpoints in close proximity are detected, all but one SVs are filtered.">\n'\
            '##FILTER=<ID=Duplicate_with_insertion,Description="Breakend SVs that are inferred to be the same as any of detected insertions">\n'\
            '##FILTER=<ID=Duplicate_with_close_insertion,Description="When multiple insertions in close proximity are detected, all but one insertions are filtered.">\n'\
            '##FILTER=<ID=SV_with_decoy,Description="SVs involving decoy contigs">\n'\
            '##FILTER=<ID=Too_small_size,Description="Insertions whose size is below the threshould (currently 100bp)">\n'\
            '##FILTER=<ID=Too_low_VAF,Description="SVs whose variant allele frequencies are inferred to be low">\n'\
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'\
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'\
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'\
            '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">\n'\
            '##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description="Length of insertion">\n'\
            '##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Sequence of insertion">\n'\
            '##ALT=<ID=DEL,Description="Deletion">\n'\
            '##ALT=<ID=INS,Description="Insertion">\n'\
            '##ALT=<ID=DUP,Description="Duplication">\n'\
            '##ALT=<ID=INV,Description="Inversion">\n'\
            '##FORMAT=<ID=TR,Number=1,Type=Integer,Description="The number of reads around the breakpoints">\n'\
            '##FORMAT=<ID=VR,Number=1,Type=Integer,Description="The number of variant supporting reads determined in the validation realignment step">'

    with open(result_file, 'r') as hin, open(output_vcf, 'w') as hout:

        dreader = csv.DictReader(hin, delimiter = '\t')
        fieldname_list = dreader.fieldnames
        is_control = True if "Checked_Read_Num_Control" in fieldname_list and "Supporting_Read_Num_Control" in fieldname_list else False

        if is_control:
            header = header + '\n' + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tCONTROL"
        else:
            header = header + '\n' + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR"

        print(header,  file = hout)

        for F in dreader:

            tchrom = F["Chr_1"]
            tid = F["SV_ID"]
            tqual = '.'
            tfilter = F["Is_Filter"]

            if F["Inserted_Seq"] != "---":
                tsvinsseq = F["Inserted_Seq"]
                tsvinslen = len(F["Inserted_Seq"])
            else:
                tsvinsseq = ''
                tsvinslen = 0
            
            tformat_sample = f'TR:VR\t{F["Checked_Read_Num_Tumor"]}:{F["Supporting_Read_Num_Tumor"]}'
            if is_control:
                tformat_sample = tformat_sample + f'\t{F["Checked_Read_Num_Control"]}:{F["Supporting_Read_Num_Control"]}' 

            if F["Chr_1"] == F["Chr_2"] and F["Dir_1"] == '+' and F["Dir_2"] == '-':

                tpos = int(F["Pos_1"])
                tref = ref_tb.fetch(tchrom, tpos - 1, tpos)
                if tref == '' or tref is None: continue
                tsvlen = int(F["Pos_2"]) - int(F["Pos_1"]) - 1
                tend = int(F["Pos_2"]) - 1
    
                # Deletion
                if tsvlen > tsvinslen:
                    talt = "<DEL>"
                    tsvlen = int(F["Pos_2"]) - int(F["Pos_1"]) - 1
                    tinfo = f"END={tend};SVTYPE=DEL;SVLEN=-{tsvlen}"

                    if tsvinslen != 0:
                        tinfo = tinfo + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                # Insertion
                elif tsvlen >= 0:
                    talt = "<INS>"
                    tinfo = f"END={tend};SVTYPE=INS;SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                else:
                    continue
                print(f"{tchrom}\t{tpos}\t{tid}\t{tref}\t{talt}\t{tqual}\t{tfilter}\t{tinfo}\t{tformat_sample}", file = hout)

            # Duplication
            elif F["Chr_1"] == F["Chr_2"] and F["Dir_1"] == '-' and F["Dir_2"] == '+' and F["Pos_1"] != '1': 

                tpos = int(F["Pos_1"])
                tref = ref_tb.fetch(tchrom, tpos - 1, tpos)
                if tref == '' or tref is None: continue
                talt = "<DUP>"
                tend = int(F["Pos_2"]) 
                tsvlen = int(F["Pos_2"]) - int(F["Pos_1"]) + 1
                tinfo = f"END={tend};SVTYPE=DUP;SVLEN={tsvlen}"
                if tsvinslen != 0:
                    tinfo = tinfo + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom}\t{tpos}\t{tid}\t{tref}\t{talt}\t{tqual}\t{tfilter}\t{tinfo}\t{tformat_sample}", file = hout)

            # Breakend 
            elif F["Chr_1"] != F["Chr_2"]:
                tchrom1 = F["Chr_1"]
                tpos1 = int(F["Pos_1"])
                tref1 = ref_tb.fetch(tchrom1, tpos1 - 1, tpos1)
                if tref1 == '' or tref1 is None: continue

                tchrom2 = F["Chr_2"]
                tpos2 = int(F["Pos_2"])
                tref2 = ref_tb.fetch(tchrom2, tpos2 - 1, tpos2)
                if tref2 == '' or tref2 is None: continue

                tbracket = ']' if F["Dir_2"] == '+' else '['
                if F["Dir_1"] == '+':
                    talt1 = f'{tref1}{tsvinsseq}{tbracket}{tchrom2}:{tpos2}{tbracket}'
                else:
                    talt1 = f'{tbracket}{tchrom2}:{tpos2}{tbracket}{tsvinsseq}{tref2}' 

                tinfo1 = f"SVTYPE=BND;MATEID={tid}_1"
                if tsvinslen != 0: tinfo1 = tinfo1 + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom1}\t{tpos1}\t{tid}_0\t{tref1}\t{talt1}\t{tqual}\t{tfilter}\t{tinfo1}\t{tformat_sample}", file = hout)

                # tchrom2 = F["Chr_2"]
                # tpos = int(F["Pos_2"])
                # tref = ref_tb.fetch(tchrom2, tpos - 1, tpos)
                # if tref == '' or tref is None: continue
                tbracket = ']' if F["Dir_1"] == '+' else '['
                tsvinsseq = reverse_complement(tsvinsseq)
                if F["Dir_2"] == '+':
                    talt2 = f'{tref2}{tsvinsseq}{tbracket}{tchrom1}:{tpos1}{tbracket}'
                else:
                    talt2 = f'{tbracket}{tchrom1}:{tpos1}{tbracket}{tsvinsseq}{tref2}'

                tinfo2 = f"SVTYPE=BND;MATEID={tid}_0"
                if tsvinslen != 0: tinfo2 = tinfo2 + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom2}\t{tpos2}\t{tid}_1\t{tref2}\t{talt2}\t{tqual}\t{tfilter}\t{tinfo2}\t{tformat_sample}", file = hout)

            else:
                tpos = int(F["Pos_1"])
                tref = ref_tb.fetch(tchrom, tpos - 1, tpos)
                if tref == '' or tref is None: continue
                talt = "<INV>"
                tend = int(F["Pos_2"]) 
                tsvlen = int(F["Pos_2"]) - int(F["Pos_1"]) + 1
                tinfo = f"END={tend};SVTYPE=INV;SVLEN={tsvlen}"
                if tsvinslen != 0:
                    tinfo = tinfo + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom}\t{tpos}\t{tid}\t{tref}\t{talt}\t{tqual}\t{tfilter}\t{tinfo}\t{tformat_sample}", file = hout)