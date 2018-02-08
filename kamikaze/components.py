from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import nt_search, seq3
from Bio.SeqRecord import SeqRecord
from .parts import PAM, Slug, EditCassette

# TODO: Finish refactoring terminology change from Payload -> Slug
class CassetteFactory():
    """
    Object for generating system components.

    Attributes
    ----------
    exposure : float
        Exposure in seconds.

    Methods
    -------
    find_pam_sites()
        Find all PAM sites
    nearest_pam_site(target)
        Find PAM site nearest to 'target'
    build_slug(target,pam_site=None)
        Build slug for target and pam_site pair.
        Will find nearest pam_site to target if not defined

    """
    def __init__(self,filename,region):
        """This is the form of a docstring
        Parameters
        ----------
        filename : str
        region : CDS region of sequence
        """

        assert isinstance(region,slice)
        self.region = region
        self.filename = filename
        self.reference = next(SeqIO.parse(filename,'fasta'))
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
        
    def build_slug(self,target,pam_site=None):
        """ Build a slug (|EDIT|...|PAM|) for a given target and pam
        uses nearest pam if not defined
        Parameters
        ----------
        filename : str
        region : CDS region of sequence
        """
        if pam_site is None:
            pam_site = self.nearest_pam_site(target)
        return Payload(self.reference.seq,target,pam_site)

    def build_edit_cassette(self,target,pam_site=None,mut='GCU'):
        pl = self.build_payload(target,pam_site)
        return EditCassette(pl,up_margin=10,down_margin=10)

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

    def gRNA_seq(self,pam_site=None):
        if pam_site is None:
            pam_site = self.nearest_pam_site()
        return self.ref_seq[self.gRNA_site(pam_site)]
    
    def build(self,donor_seq,pam_site=None):
        if pam_site is None:
            pam_site = self.nearest_pam_site()
        grna_seq = self.gRNA_seq(pam_site)
        pam_seq = self.ref_seq[pam_site]
        repair = self.repair_seq(donor_seq,pam_site)
        return Seq(str(grna_seq)+str(pam_seq)+str(repair),IUPAC.unambiguous_dna)
