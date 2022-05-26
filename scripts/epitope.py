from pyensembl import EnsemblRelease
from common import Variant
from Bio.Seq import translate

def translate_dna(seq):
    return translate(seq, to_stop=True)

def create_epitope(variants):
    ens = EnsemblRelease(106)
    for i, variant in enumerate(variants):
        if ':n.' in str(variant.hgvsc):
            del variants[i]
        variant_start = variant.hgvsc.posedit.pos.start.base - 1
        variant_end = variant.hgvsc.posedit.pos.end.base - 1
        transcript = variant.transcript.split('.')[0]
        transcript_ensembl = ens.transcript_by_id(transcript)
        try:
            sequence = transcript_ensembl.coding_sequence
        except ValueError:
            sequence = transcript_ensembl.sequence
        if variant.alt == 'DEL':
            aa_pos = round(variant_start/3)
            start = aa_pos - 13
            cDNA = sequence[:variant_start] + sequence[variant_end:]
            wt_aa = translate_dna(sequence)
            mut_aa = translate_dna(cDNA)
            variant.wt = wt_aa[start:]
            variant.mut = mut_aa[start:]
        elif variant.alt == 'INS':
            continue