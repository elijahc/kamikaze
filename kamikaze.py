import pybedtools as pbt
from pybedtools import BedTool
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils.CodonUsage import SynonymousCodons
from Bio.SeqUtils import nt_search, seq3

class PAM():
    def __init__(self,strain=None, seq=None):
        if strain is not None and seq is None:
            self.strain = strain
            if strain.lower()=='cas9':
                self.seq = Seq('NGG',IUPAC.unambiguous_dna)

        elif strain is None and seq is not None:
            self.seq = seq
            if self.seq.upper()=='NGG': 
                self.strain='cas9'

    def find_in(self,search_seq):
        return nt_search(str(search_seq),str(self.seq))[1:]

def codon_synonyms(codon):
    AAA = seq3(codon.translate()).upper()
    return SynonymousCodons[AAA]

class EditingCassette():
    def __init__(self,ref_seq,target,pam,up_margin=10,down_margin=10):
        self.ref_seq = ref_seq
        self.up_margin = up_margin
        self.down_margin = down_margin

        if isinstance(target,int):
            self.target = slice(target,target+3)
        elif isinstance(target,tuple):
            self.target = slice(target[0],target[1])
        elif isinstance(target,slice):
            self.target = target
        else: raise ValueError

        self.pam = pam
        self.pam_sites = [slice(s,s+3) for s in self.pam.find_in(self.ref_seq)]

    def nearest_pam_site(self):
        site_stop = min([s.stop for s in self.pam_sites if s.stop<self.target.start],key=lambda x:abs(x-self.target.start))
        return slice(site_stop-3,site_stop)

    def upstream_HA(self,pam_site):
        uha_start = max([0,pam_site.start-self.up_margin])
        return self.ref_seq[uha_start:pam_site.start]

    def downstream_HA(self,pam_site):
        dha_stop = min([len(self.ref_seq),self.target.stop+self.down_margin])
        return self.ref_seq[self.target.stop:dha_stop]

    def gRNA_site(self,pam_site=None):
        if pam_site is None:
            pam_site = self.nearest_pam_site()
        return slice(max([pam_site.start-20,0]),pam_site.start)

    def gRNA_seq(self,pam_site=None):
        if pam_site is None:
            pam_site = self.nearest_pam_site()
        return self.ref_seq[self.gRNA_site(pam_site)]
    
    def repair_seq(self,donor_seq,pam_site):
        pam_seq = self.ref_seq[pam_site]
        syn_pams = codon_synonyms(pam_seq)
        uha = self.upstream_HA(pam_site)
        dha = self.downstream_HA(pam_site)
        edit_window = str(syn_pams[0])+str(self.ref_seq[pam_site.stop:self.target.start])+str(donor_seq)
        out_seq = str(uha).lower()+str(edit_window).upper()+str(dha).lower()
        return Seq(out_seq,IUPAC.unambiguous_dna)

    def build(self,donor_seq,pam_site=None):
        if pam_site is None:
            pam_site = self.nearest_pam_site()
        grna_seq = self.gRNA_seq(pam_site)
        pam_seq = self.ref_seq[pam_site]
        repair = self.repair_seq(donor_seq,pam_site)
        return Seq(str(grna_seq)+str(pam_seq)+str(repair),IUPAC.unambiguous_dna)

def gen_codon_bed(cds_start,cds_stop,chrom='chr1'):
    idxs = np.arange(cds_start,cds_stop,3)
    intervals = [ pbt.create_interval_from_list([chrom,str(i),str(i+3)]) for i in idxs]
    return pbt.BedTool(i for i in intervals)

def load_cds(fi='kras_cds.fa'):
    return next(SeqIO.parse(fi,'fasta'))

def gen_target_HA(dna_seq,bp_idx,margin=5):

    return dna_seq[bp_idx-margin:bp_idx+margin]