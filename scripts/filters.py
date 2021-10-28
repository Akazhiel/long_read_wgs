"""
@author: Jonatan González Rodríguez <jonatan.gonzalez.r@outlook.com>
"""

import vcfpy


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
        if tumor_DP >= 10 and tumor_DV >= 7 and normal_DV <= 1:
            writer.write_record(record)


def filter_callers(inp, out, num_callers):
    reader = vcfpy.Reader.from_path(inp)
    writer = vcfpy.Writer.from_path(out, reader.header)
    for record in reader:
        called = [x for x in record.calls if x.data['DR'] is not None and 'Normal' not in x.sample]
        if len(called) >= num_callers:
            writer.write_record(record)
