"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
"""

import pandas as pd
import vcfpy


def filter_somatic(inp, out, caller):
    reader = vcfpy.Reader.from_path(inp)
    writer = vcfpy.Writer.from_path(out, reader.header)
    for record in reader:
        calls = {x.sample: x.data for x in record.calls}
        token = calls[f'{caller}_Normal']['AD']
        normal_DR = token if type(token) is int else token[0]
        token = calls[f'{caller}_Tumor']['AD']
        tumor_DR = token if type(token) is int else token[0]
        token = calls[f'{caller}_Normal']['AD']
        normal_DV = token if type(token) is int else token[1]
        token = calls[f'{caller}_Tumor']['AD']
        tumor_DV = token if type(token) is int else token[1]
        tumor_DP = tumor_DR + tumor_DV
        normal_DP = normal_DR + normal_DV
        tumor_VAF = round((tumor_DV/tumor_DP) * 100, 3) if tumor_DP != 0 else 0
        normal_VAF = round((normal_DV/normal_DP) * 100, 3) if normal_DP != 0 else 0
        t2n_ratio = tumor_VAF/normal_VAF if normal_VAF != 0 else 5
        if tumor_DP >= 10 and tumor_DV >= 7 and t2n_ratio >= 5 and normal_DV <= 3:
            writer.write_record(record)


def filter_callers(inp, out, num_callers):
    reader = vcfpy.Reader.from_path(inp)
    writer = vcfpy.Writer.from_path(out, reader.header)
    for record in reader:
        called = [x for x in record.calls if x.data['AD'] is not None and 'Normal' not in x.sample]
        if len(called) >= num_callers:
            writer.write_record(record) 

def prioritize_variants(annot_sv_df):
    variants = pd.read_csv(annot_sv_df, sep = "\t").iloc[:, :35]
    variants[["Location1", "Location2"]] = variants["Location"].str.split('-', expand=True)
    variants = variants[variants.Location != 'txStart-txEnd']
    variants = variants[~variants.apply(lambda x: x['Location1'] == x['Location2'] and 'intron' in x['Location'], axis = 1)]
    variants = variants[~variants.apply(lambda x: x['SV_type'] == 'DEL' and 'txStart' in x['Location1'] , axis = 1)]
    variants['Priority'] = variants.apply(lambda x: 2 if (x['SV_type'] == 'DEL' and 'txEnd' in x['Location2']) else 3 if (x['SV_type'] not in ['INS', 'DEL']) else 1, axis = 1)
    variants['SV_chrom'] = 'chr' + variants['SV_chrom'].astype(str)
    variants['Frameshift'] = variants['Frameshift'].map({"yes":True,"no":False})
    variants.to_csv('annotsv_prioritized.tsv', sep="\t", index=False)
    return variants
