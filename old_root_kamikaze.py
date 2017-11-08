import pybedtools as pbt
from pybedtools import BedTool
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.CodonUsage import SynonymousCodons
from Bio.SeqUtils import nt_search, seq3

def codon_synonyms(codon):
    AAA = seq3(codon.translate()).upper()
    return SynonymousCodons[AAA]

class Variant():
    def __init__(self,ref_seq,target,insert,pam_site):
        self.target = target
        self.insert = insert
        self.pam_site = pam_site

def load_cds(fi='kras_cds.fa'):
    return next(SeqIO.parse(fi,'fasta'))
