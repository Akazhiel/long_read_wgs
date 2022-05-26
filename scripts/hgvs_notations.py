# import hgvs.assemblymapper
# import hgvs.parser
# import pandas as pd
# from cdot.hgvs.dataproviders import JSONDataProvider
# from hgvs.transcriptmapper import HGVSUsageError
# from scripts.common import Variant

# # Load chromosome map to map UCSC identifiers to RefSeq

# # Have to add internal path for docker

# chrmap = pd.read_csv("./GRCh38_RefSeq2UCSC.txt", sep = "\t", names=['RefSeq', 'UCSC'], index_col=1).to_dict()['RefSeq']

# # Include path inside the docker where this will be located.

# local_json = ["cdot-0.2.6.ensembl.grch38.json.gz", "cdot-0.2.6.refseq.grch38.json.gz"]

# # Initialize parser and variant mapper.

# parse = hgvs.parser.Parser()
# hdp = JSONDataProvider(local_json)
# vm = hgvs.assemblymapper.AssemblyMapper(
#     hdp, assembly_name='GRCh38', alt_aln_method='splign', replace_reference=True)

# # Add HGVSg notation for variants, insertions are marked with an A to circumvent error when using N base

# def add_variant_hgvsg(variants_df):
#     variants = list()
#     for _, x in variants_df.iterrows():
#         variant            = Variant()
#         variant.chrom      = x['SV_chrom']
#         variant.start      = x['SV_start']
#         variant.end        = x['SV_end']
#         variant.alt        = x['SV_type']
#         variant.gene       = x['Gene_name']
#         variant.transcript = x['Tx']
#         variant.priority   = x['Priority']
#         variant.frameshift = x['Frameshift']
#         if 'INS' in variant.alt:
#             variant.hgvsg = parse.parse_hgvs_variant(chrmap[x['SV_chrom']] + ':g.' + str(x['SV_start']) + '_' +  str(x['SV_end']) + 'insA')
#         else:
#             variant.hgvsg = parse.parse_hgvs_variant(chrmap[x['SV_chrom']] + ':g.' + str(x['SV_start']) + '_' +  str(x['SV_end']) + 'del')
#         variants.append(variant)

# def add_variant_hgvsc(variants):
#     for variant in variants:
#         try:
#             hgvs_cds = vm.g_to_c(variant.hgvsg, variant.transcript)
#             variant.hgvsc = hgvs_cds
#         except HGVSUsageError:
#             hgvs_cds = vm.g_to_n(variant.hgvsg, variant.transcript)
#             variant.hgvsc = hgvs_cds
#     return variants

# def add_variant_hgvs(annotsv_df):
#     variant_hgvsg = add_variant_hgvsg(annotsv_df)
#     variant_hgvsc = add_variant_hgvsc(variant_hgvsg)
#     return variant_hgvsc