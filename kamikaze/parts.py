from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

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

class Payload():
    def __init__(self,ref_seq,target,pam_site):
        assert isinstance(pam_site,SeqFeature)
        self.ref_seq = ref_seq
        self.target = target
        self.pam_site = pam_site
    
    def __str__(self):
        start = min([self.target.start,self.pam_site.location.start])
        stop = max([self.target.stop,self.pam_site.location.end])
        return str(self.ref_seq[start:stop])