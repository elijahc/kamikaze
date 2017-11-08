from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import nt_search, seq3
from Bio.SeqRecord import SeqRecord
from .parts import PAM, Payload

class CassetteFactory():
    def __init__(self,fi,region):
        assert isinstance(region,slice)
        self.region = region
        self.fi = fi
        self.reference = next(SeqIO.parse(fi,'fasta'))
        self.reference.features = [
            SeqFeature(FeatureLocation(region.start,region.stop),type='cds',strand=1)
            ]
        self.pam = PAM(strain='cas9')
        self._feature_by_type = lambda type: [f for f in self.reference.features if f.type==type]

    def find_pam_sites(self):
        search_space = self.reference.seq
        pam_seq = self.pam.seq
        pam_len = len(pam_seq)
        fwd_hits = nt_search(str(search_space),str(pam_seq))[1:]
        rev_hits = nt_search(str(search_space.complement()),str(pam_seq)[::-1])[1:]
        i2ps = lambda i,s: SeqFeature(FeatureLocation(i,i+len(pam_seq)),type='pam_site',strand=s)
        pam_sites = [i2ps(i,1) for i in fwd_hits]
        pam_sites.extend([i2ps(i,-1) for i in rev_hits])
        # pam_recs = SeqRecord(search_space,name='pam sites',features=pam_sites)
        return pam_sites
    
    def nearest_pam_site(self,target):
        dist = lambda t,ps: min([abs(t.start-ps.location.end), abs(t.stop-ps.location.start)])
        exclusive = lambda t,ps: ps.location.end < t.start or ps.location.start > t.stop
        pam_sites = [ps for ps in self.find_pam_sites() if exclusive(target,ps)]
        pam_distances = [dist(target,ps) for ps in pam_sites]
        return min(zip(pam_sites,pam_distances),key=lambda x:x[1])[0]
        
    def build_payload(self,target,pam_site=None):
        if pam_site is None:
            pam_site = self.nearest_pam_site(target)
        return Payload(self.reference.seq,target,pam_site)

    def build_edit_cassette(self,target,pam_site=None,mut='GCU'):
        pl = self.build_payload(target,pam_site)
        ref_seq = self.reference.seq
        edit_pl = ''
        if pl.target.start < pl.pam_site.location.start:
            edit_pl+=mut+str(ref_seq[pl.target.end:pl.pam_site.location.end]
        else:
            edit_pl+=str(ref_seq[pl.pam_site.location.start:pl.target.start])+mut

        return edit_pl


# TODO Rewrite Editing Cassette to support CassetteFactory
class EditingCassette():
    def __init__(self,target,pam,up_margin=10,down_margin=10):
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