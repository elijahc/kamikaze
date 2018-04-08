from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from .edits import CodonEdit,Swap
from .pams import *

class Slug():
    """Target codon, pam, and any gap between

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
        self.target_seq = self.ref_seq[self.target]
        self.pam_site = pam_site
        pam_slc = slice(self.pam_site.location.start,self.pam_site.location.end)
        self.pam_seq = self.ref_seq[pam_slc]
        self.start = min([self.target.start,self.pam_site.location.start])
        self.stop = max([self.target.stop,self.pam_site.location.end])
        
        if self.target.start < self.pam_site.location.start:
            # [Target]-[gap]-[PAM]
            self.pam_rel_loc = 'downstream'
        else:
            # [PAM]-[gap]-[Target]
            self.pam_rel_loc = 'upstream'

    def print(self):
        out = ''
        gap = ['-']*(len(self)-6)
        if self.pam_rel_loc == 'downstream':
            out+='TGT'
            out+=''.join(gap)
            out+='PAM'
        else:
            out+='PAM'
            out+=''.join(gap)
            out+='TGT'
        print(out)
        print(str(self))

    
    def __str__(self):
        return str(self.ref_seq[self.start:self.stop])

    def __len__(self):
        return len(self.ref_seq[self.start:self.stop])


class Payload():
    def __init__(self,slug,up_margin=10,down_margin=10):
        if not isinstance(up_margin,int):
            raise TypeError('up_margin and down_margin must be ints')
        elif not isinstance(down_margin,int):
            raise TypeError('up_margin and down_margin must be ints')

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
            edit_pl = mut()+str(ref_seq[pl.target.stop:pl.pam_site.location.end])
        else:
            # [PAM]-[gap]-[MUT]
            edit_pl = str(ref_seq[pl.pam_site.location.start:pl.target.start])+mut()
        
        uha = self.upstream_HA()
        dha = self.downstream_HA()

        out_seq = str(uha).upper()+str(edit_pl)+str(dha).upper()
        return Seq(out_seq,IUPAC.unambiguous_dna)

class EditCassette():
    def __init__(self,slug,sg_promoter,
                 chip_primer='GGGTTTGAAGGATACCAGCT',
                 up_margin=10,
                 down_margin=10,
                 crispr_len=20,
                 name='<default>',
                 description='<default>'):
        self.slug = slug
        self.name = name
        self.description = description
        self.chip_primer = chip_primer
        self.sg_promoter = sg_promoter
        self.up_margin = up_margin
        self.down_margin = down_margin
        self.crispr_len = crispr_len

        if isinstance(sg_promoter,Seq):
            pass
        elif isinstance(sg_promoter,str):
            gRNA_promoter = Seq(sg_promoter,IUPAC.unambiguous_dna)
        else:
            raise TypeError("sg_promoter must be a string or Seq object")

    def payload(self):
        return Payload(self.slug,self.up_margin,self.down_margin)

    def gRNA_site(self):
        """ Location of gRNA in refseq
        ...
        Returns
        ----------
        slice
            slice location of gRNA_site within slug.ref_seq
        """
        
        pam_site_loc = self.slug.pam_site.location
        if pam_site_loc.strand == 1:
            return slice(max([pam_site_loc.start-self.crispr_len,0]),int(pam_site_loc.start))
        else:
            return slice(int(pam_site_loc.end),min([len(self.slug.ref_seq),pam_site_loc.end+self.crispr_len]))
    
    def gRNA_seq(self):
        pam_site_loc = self.slug.pam_site.location
        seq = self.slug.ref_seq[self.gRNA_site()]
        if pam_site_loc.strand == -1: 
            seq = seq.complement()

        return seq

    def assemble_oligo(self,sp_primer,edit=None):
        """ Assemble completed edit cassette oligo
        ...
        Parameters
        ----------
        sp_primer : Seq or str
            Subpool Primer sequence (~20bp)

        gRNA_promoter : Seq or str
            Promoter sequence for the CRISPR guide RNA (~35bp)

        edit : str or CodonEdit-like
            if a string is passed it will substitute the codon for the string passed

        Returns
        ----------
        SeqRecord
            SeqRecord with all relevant features annotated
        """
        if isinstance(sp_primer,Seq):
            pass
        elif isinstance(sp_primer,str):
            sp_primer = Seq(sp_primer,IUPAC.unambiguous_dna)
        else:
            raise TypeError("sp_primer must be a string or Seq object")

        if isinstance(edit,CodonEdit):
            mut = edit
        elif isinstance(edit,str):
            tgt_slc = self.slug.target
            tgt_codon = self.slug.ref_seq[tgt_slc]
            mut = Swap(codon=tgt_codon,new_codon=edit)
        else:
            raise TypeError("edit must be of type str or CodonEdit")

        slug_len = len(self.slug)
        pl = self.payload().assemble(mut)
        self.crispr_rna = self.gRNA_seq()

        parts = [sp_primer,pl,self.sg_promoter,self.crispr_rna,self.chip_primer]
        labels = ['subpool_primer','payload','gRNA_promoter','crRNA','end-primer']

        features = []
        start=0
        for lab,part in zip(labels,parts):
            loc = FeatureLocation(start,start+len(part))
            sf = SeqFeature(type=lab,location=loc)
            features.append(sf)
            start +=len(part)

        out_seq = sp_primer+pl+self.sg_promoter+self.crispr_rna+self.chip_primer
        annotations = {
            'target':self.slug.target,
            'slug_len':len(self.slug)
        }

        desc = "|tgt:"+str(self.slug.target.start)+".."+str(self.slug.target.stop)+"|"+str(mut)
        edit_seq_rec = SeqRecord(seq=out_seq,
                                 id=self.name+'.'+str(id(out_seq)),
                                 name=self.name,
                                 description=desc,
                                 annotations=annotations,
                                 features=features)
        return edit_seq_rec
