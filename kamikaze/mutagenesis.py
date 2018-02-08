from .components import CassetteFactory

class Mutagenesis(CassetteFactory):
    """Base class for generating cassettes codon by codon

    ...

    Attributes
    ----------
    region  : slice
        Slice object representing CDS of sequence
    targets : list
        List of Slice objects for each codon

    Methods
    -------
    gen_slugs()
        Returns a list of slugs for each target in self.targets
    gen_payloads()
        Returns a list of payloads for each target in self.targets
    gen_cassettes()
        Returns a list of edit cassettes for each target in self.targets
    """
    def __init__(self,filename,region,**kwargs):
        super().__init__(filename,region,**kwargs)
        targets = []
        self.slugs = None
        self.payloads = None
        self.cassettes = None

        i = self.region.start
        while i+3 < self.region.stop:
            targets.append(slice(i,i+3))
            i+=3
        self.targets = targets 
    
    def gen_slugs(self):
        self.slugs = [self.build_slug(tgt) for tgt in self.targets]
        return self.slugs

    def gen_payloads(self):
        if self.slugs is not None:
            payloads = [self.build_payload(slug) for slugs in self.slugs]
        else:
            payloads = []
            for tgt in self.targets:
                slug = self.build_slug(tgt)
                pl = self.build_payload(slug)
                payloads.append(pl)
            
        self.payloads = payloads

        return self.payloads

    def gen_cassettes(self,sg_promoter,**kwargs):
        if self.slugs is not None:
            cassettes = [self.build_cassette(slug,**kwargs) for slugs in self.slugs]
        else:
            cassettes = []
            for tgt in self.targets:
                slug = self.build_slug(tgt)
                cass = self.build_cassette(slug,sg_promoter,**kwargs)
                cassettes.append(cass)
            
        self.cassettes = cassettes

        return self.cassettes

class AlanineScan(Mutagenesis):
    """CassetteFactory for generating Alanine Scans

    ...

    Attributes
    ----------
    region  : slice
        Slice object representing CDS of sequence
    targets : list
        List of Slice objects for each codon

    Methods
    -------
    gen_slugs()
        Returns a list of slugs for each target in self.targets
    gen_payloads()
        Returns a list of payloads for each target in self.targets
    gen_cassettes()
        Returns a list of edit cassettes for each target in self.targets
    """
    def __init__(self,filename,region,**kwargs):
        super().__init__(filename,region,**kwargs)
    
    def assemble_oligos(self,sp_primer,sg_promoter,crispr_rna_len=20,max_len=230):
        oligos = []
        if self.cassettes is not None:
            for ec in self.cassettes:
                oligo = ec.assemble_oligo(sp_primer=sp_primer,edit='GCU')
                oligos.append(oligo)
        else:
            for tgt in self.targets:
                slug = self.build_slug(tgt)

                if self.ha_margins=='max':
                    # spp + HA/2 + slug_len + HA/2 + sg_promoter + crRNA + ... = max_len
                    len_static = len(sp_primer) + len(slug) +len(sg_promoter) + crispr_rna_len + 20
                    ha_margins = (max_len - len_static)/2
                else:
                    ha_margins=self.ha_margins

                self.cassette_args = dict(
                    slug=slug,
                    sg_promoter=sg_promoter,
                    up_margin=int(ha_margins),
                    down_margin=int(ha_margins),
                    crispr_len=crispr_rna_len
                )

                ec = self.build_cassette(**self.cassette_args)

                oligo = ec.assemble_oligo(sp_primer=sp_primer,edit='XXX')
                oligos.append(oligo)

        return oligos
    
# TODO: class CodonSubstitution()
# TODO: class TruncationScan()
# TODO: class SingleBaseDeletion()
# TODO: class DoubleBaseDeletion()
# TODO: class CodonDeletionScan()