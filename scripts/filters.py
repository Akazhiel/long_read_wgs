"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
"""

import vcfpy
import numpy as np

def filter_somatic(inp, out, caller):
    reader = vcfpy.Reader.from_path(inp)
    writer = vcfpy.Writer.from_path(out, reader.header)
    for record in reader:
        calls = {x.sample: x.data for x in record.calls}
        token = calls[f'{caller}_Normal']['DR']
        normal_DR = token if type(token) is int else 0
        token = calls[f'{caller}_Tumor']['DR']
        tumor_DR = token if type(token) is int else 0
        token = calls[f'{caller}_Normal']['DV']
        normal_DV = token if type(token) is int else 0
        token = calls[f'{caller}_Tumor']['DV']
        tumor_DV = token if type(token) is int else 0
        tumor_DP = tumor_DR + tumor_DV
        normal_DP = normal_DR + normal_DV
        tumor_VAF = np.around(tumor_DV/tumor_DP, 3) if tumor_DP != 0 else 0
        normal_VAF = np.around(normal_DV/normal_DP, 3) if normal_DP != 0 else 0
        t2n_ratio = tumor_VAF/normal_VAF if normal_VAF != 0 else 5
        if tumor_DP >= 10 and tumor_DV >= 7 and t2n_ratio >= 5:
            writer.write_record(record)


def filter_callers(inp, out, num_callers):
    reader = vcfpy.Reader.from_path(inp)
    writer = vcfpy.Writer.from_path(out, reader.header)
    for record in reader:
        called = [x for x in record.calls if x.data['DR'] is not None and 'Normal' not in x.sample]
        if len(called) >= num_callers:
            writer.write_record(record)
