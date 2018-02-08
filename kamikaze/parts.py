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

class Slug():
    """
    Target codon, pam, and any gap between

    Attributes
    ----------
    ref_seq : reference
        
    target : slice
        
    pam_site : SeqFeature

    start : int

    stop : int

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
    def __init__(self,ref_seq,target,pam_site):
        assert isinstance(pam_site,SeqFeature)
        self.ref_seq = ref_seq
        self.target = target
        self.pam_site = pam_site
        self.start = min([self.target.start,self.pam_site.location.start])
        self.stop = max([self.target.stop,self.pam_site.location.end])
        
        if self.target.start < self.pam_site.location.start:
            # [Target]-[gap]-[PAM]
            self.pam_rel_loc = 'downstream'
        else:
            # [PAM]-[gap]-[Target]
            self.pam_rel_loc = 'upstream'
    
    def __str__(self):
        return str(self.ref_seq[self.start:self.stop])

class Payload():
    def __init__(self,slug,up_margin=10,down_margin=10):
        self.slug = slug
        self.up_margin = up_margin
        self.down_margin = down_margin

    def upstream_HA(self):
        uha_start = max([0,self.slug.start-self.up_margin])
        return self.slug.ref_seq[uha_start:self.slug.start]

    def downstream_HA(self):
        dha_stop = min([len(self.slug.ref_seq),self.slug.stop+self.down_margin])
        return self.slug.ref_seq[self.slug.stop:dha_stop]

    def assemble(self,mut):
        # pam_seq = self.payload.ref_seq[pam_site]
        # syn_pams = codon_synonyms(pam_seq)
        pl = self.slug
        ref_seq = pl.ref_seq
        edit_payload = ''
        if pl.pam_rel_loc == 'downstream':
            # [MUT]-[gap]-[PAM]
            edit_pl = mut+str(ref_seq[pl.target.stop:pl.pam_site.location.end])
        else:
            # [PAM]-[gap]-[MUT]
            edit_pl = str(ref_seq[pl.pam_site.location.start:pl.target.start])+mut
        
        uha = self.upstream_HA()
        dha = self.downstream_HA()

        out_seq = str(uha).upper()+str(edit_pl)+str(dha).upper()
        return Seq(out_seq,IUPAC.unambiguous_dna)

class EditCassette():
    def __init__(self,slug,up_margin=10,down_margin=10,grna_len=20):
        self.slug = slug
        self.up_margin = up_margin
        self.down_margin = down_margin
        self.grna_len = grna_len

    def payload(self):
        return Payload(slug,self.up_margin,self.down_margin)

    def gRNA_site(self):
        pam_site_loc = self.payload.pam_site.location
        if pam_site_loc.strand == 1:
            return slice(max([pam_site_loc.start-grna_len,0]),int(pam_site_loc.start))
        else:
            return slice(int(pam_site_loc.start),min([len(self.payload.ref_seq),pam_site_loc.end+grna_len]))


    # def __str__(self):
    #     return self.assemble_oligo(mut='GCU')