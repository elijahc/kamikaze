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
        self.start = min([self.target.start,self.pam_site.location.start])
        self.stop = max([self.target.stop,self.pam_site.location.end])
    
    def __str__(self):
        return str(self.ref_seq[self.start:self.stop])

class EditCassette():
    def __init__(self,payload,up_margin=10,down_margin=10):
        self.payload = payload
        self.up_margin = up_margin
        self.down_margin = down_margin

    def gRNA_site(self,pam_site=None):
        if pam_site is None:
            pam_site = self.nearest_pam_site()
        return slice(max([pam_site.start-20,0]),pam_site.start)

    def upstream_HA(self):
        uha_start = max([0,self.payload.start-self.up_margin])
        return self.payload.ref_seq[uha_start:self.payload.start]

    def downstream_HA(self):
        dha_stop = min([len(self.payload.ref_seq),self.payload.stop+self.down_margin])
        return self.payload.ref_seq[self.payload.stop:dha_stop]

    def repair_payload(self,mut):
        # pam_seq = self.payload.ref_seq[pam_site]
        # syn_pams = codon_synonyms(pam_seq)
        pl = self.payload
        ref_seq = pl.ref_seq
        edit_payload = ''
        if pl.target.start < pl.pam_site.location.start:
            edit_pl = mut+str(ref_seq[pl.target.stop:pl.pam_site.location.end])
        else:
            edit_pl = str(ref_seq[pl.pam_site.location.start:pl.target.start])+mut

        return edit_pl

    def repair_seq(self,mut='GCU'):
        edited_payload = self.repair_payload(mut)
        uha = self.upstream_HA()
        dha = self.downstream_HA()
        # edit_window = str(syn_pams[0])+str(self.ref_seq[pam_site.stop:self.target.start])+str(donor_seq)
        out_seq = str(uha).upper()+str(repair_seq).upper()+str(dha).upper()
        return Seq(out_seq,IUPAC.unambiguous_dna)

    # def __str__(self):
    #     return self.assemble_oligo(mut='GCU')