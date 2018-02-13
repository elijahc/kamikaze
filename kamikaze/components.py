from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import nt_search, seq3
from Bio.SeqRecord import SeqRecord
from .parts import PAM, Slug, Payload, EditCassette

# TODO: Finish refactoring terminology change from Payload -> Slug
class CassetteFactory():
    """
    Object for generating system components.

    Attributes
    ----------
    region : slice
        Slice of the mRNA cds
    filename : string
        path sequence fasta
    reference : SeqIO
    pam : PAM
    ha_margins : int

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
    def __init__(self,filename,region,lib_name='default',ha_margins=10):

        assert isinstance(region,slice)
        self.region = region
        self.filename = filename
        self.pam_sites = None
        self.reference = next(SeqIO.parse(filename,'fasta'))
        self.ha_margins = ha_margins
        self.reference.features = [
            SeqFeature(FeatureLocation(region.start,region.stop),type='cds',strand=1)
            ]
        self.pam = PAM(strain='cas9')
        self.lib_name = lib_name
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
        self.pam_sites = pam_sites
        return pam_sites
    
    def nearest_pam_site(self,target,n=1):
        """ Find nearest pam_site to target

        Parameters
        ----------
        target : slice
        """

        # Shortcut for calculating distance between 2 regions (target and pam_site?)
        dist = lambda t,ps: min([abs(t.start-ps.location.end), abs(t.stop-ps.location.start)])
        exclusive = lambda t,ps: ps.location.end < t.start or ps.location.start > t.stop
        if self.pam_sites is None:
            pam_sites = [ps for ps in self.find_pam_sites() if exclusive(target,ps)]
        else:
            pam_sites = [ps for ps in self.pam_sites if exclusive(target,ps)]
        pam_distances = [dist(target,ps) for ps in pam_sites]
        sorted_pams = sorted(zip(pam_sites,pam_distances),key=lambda x:x[1])
        return [p[0] for p in sorted_pams[:n]]

        
    def build_slug(self,target,pam_site=None,n=1,max_len=25):
        """ Build a slug (|EDIT|...|PAM|) for a given target and nearest pam or specified one
        ...
        Parameters
        ----------
        target : slice
        pam_site : SeqFeature
        """
        if pam_site is None:
            pam_sites = self.nearest_pam_site(target,n)
        slugs = [Slug(self.reference.seq,target,ps) for ps in pam_sites]
        return [s for s in slugs if len(s)<max_len]

    def build_payload(self,slug,mut='GCU'):
        pl = Payload(slug)
        return pl

    def build_cassette(self,slug,sg_promoter,
    up_margin=10,down_margin=10,crispr_len=20,name='cassette',description='cassette'):
        """ Build a cassette
            ...
            Parameters
            ----------
            slug : Slug
            sg_promoter : str or Seq
            up_margin : int, optional
            down_margin : int, optional
            crispr_len : int, optional

            Returns
            ----------
            EditCassette
                EditCassette that still needs assembled
        """
        ec = EditCassette(slug,sg_promoter=sg_promoter,
                          up_margin=up_margin,
                          down_margin=down_margin,
                          crispr_len=crispr_len,
                          name=name,
                          description=description)
        return ec