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
        if ((tumor_DP and tumor_DV) >= 1 and normal_DP == 0) or (
            (normal_DR and tumor_DR and tumor_DV) >= 1 and normal_DV <= 1
        ):
            writer.write_record(record)
