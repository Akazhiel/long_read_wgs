import hgvs.assemblymapper
import hgvs.parser
import pandas as pd
from cdot.hgvs.dataproviders import JSONDataProvider
from hgvs.transcriptmapper import HGVSUsageError
from scripts.common import Variant

# Load chromosome map to map UCSC identifiers to RefSeq

# Have to add internal path for docker

chrmap = pd.read_csv("./GRCh38_RefSeq2UCSC.txt", sep = "\t", names=['RefSeq', 'UCSC'], index_col=1).to_dict()['RefSeq']

# Include path inside the docker where this will be located.

local_json = ["cdot-0.2.11.ensembl.grch38.json.gz", "cdot-0.2.11.refseq.grch38.json.gz"]

# Initialize parser and variant mapper.

parse = hgvs.parser.Parser()
hdp = JSONDataProvider(local_json)
vm = hgvs.assemblymapper.AssemblyMapper(
    hdp, assembly_name='GRCh38', alt_aln_method='splign', replace_reference=False, normalize=False)

# Add HGVSg notation for variants, insertions are marked with an A to circumvent error when using N base

def add_variant_hgvs(variants_df):
    variants = list()
    for _, x in variants_df.iterrows():
        variant            = Variant()
        variant.chrom      = x['SV_chrom']
        variant.start      = x['SV_start']
        variant.end        = x['SV_end']
        variant.alt        = x['SV_type']
        variant.gene       = x['Gene_name']
        variant.transcript = x['Tx']
        variant.priority   = x['Priority']
        variant.frameshift = None
        variant.length     = x['SV_length']
        variant.location   = x['Location']
        if 'INS' in variant.alt:
            variant.hgvsg = parse.parse_hgvs_variant(chrmap[x['SV_chrom']] + ':g.' + str(x['SV_start']) + '_' +  str(x['SV_end']) + 'insA')
            variant.hgvsc = add_variant_hgvsc(variant.hgvsg, variant.transcript)
        else:
            variant.hgvsg = parse.parse_hgvs_variant(chrmap[x['SV_chrom']] + ':g.' + str(x['SV_start']) + '_' +  str(x['SV_end']) + 'del')
            variant.hgvsc = add_variant_hgvsc(variant.hgvsg, variant.transcript)
        variants.append(variant)

def add_variant_hgvsc(variant_hgvsg, transcript):
    try:
        variant_hgvsc = vm.g_to_c(variant_hgvsg, transcript)
    except HGVSUsageError:
        variant_hgvsc = vm.g_to_n(variant_hgvsg, transcript)
    return variant_hgvsc