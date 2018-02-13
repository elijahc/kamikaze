from .components import CassetteFactory
from Bio.Data.CodonTable import standard_dna_table,unambiguous_dna_by_name

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
        slugs = []
        for tgt in self.targets:
            slugs.extend(self.build_slug(tgt,n=3))
        self.slugs = slugs
        return self.slugs

    def gen_payloads(self):
        if self.slugs is not None:
            payloads = [self.build_payload(slug) for slugs in self.slugs]
        else:
            payloads = []
            for tgt in self.targets:
                slugs = self.build_slug(tgt)
                pls = [self.build_payload(slug) for slug in slugs]
                payloads.extend(pls)
            
        self.payloads = payloads

        return self.payloads

    def gen_cassettes(self,sg_promoter=None,**kwargs):
        if sg_promoter is None:
            print('No guide promoter specified, using default')
            self.sg_promoter = 'TTGACAGCTAGCTCAGTCCTAGGTATAATACTAGT'
        if self.slugs is not None:
            cassettes = [self.build_cassette(slug,**kwargs) for slugs in self.slugs]
        else:
            cassettes = []
            for tgt in self.targets:
                tgt_slugs = self.build_slug(tgt,n=4)
                for slug in tgt_slugs:
                    cass = self.build_cassette(slug,self.sg_promoter,**kwargs)
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
    def __init__(self,filename,region,lib_name='Ala_scan',codon_table='Standard',**kwargs):
        super().__init__(filename,region,**kwargs)
        if codon_table not in list(unambiguous_dna_by_name.keys()):
            print(codon_table,' not in dataset')
            print('Alternative for Standard:')
            print(unambiguous_dna_by_name.keys())
            raise
        self.ala_codon = unambiguous_dna_by_name[codon_table].back_table['A']
    
    def assemble_oligos(self,sp_primer,sg_promoter=None,crispr_rna_len=20,max_len=227):
        if sg_promoter is None:
            print('No guide promoter specified, using default')
            self.sg_promoter = 'TTGACAGCTAGCTCAGTCCTAGGTATAATACTAGT'
        else:
            self.sg_promoter = sg_promoter
        oligos = []
        if self.cassettes is not None:
            for ec in self.cassettes:
                oligo = ec.assemble_oligo(sp_primer=sp_primer,sg_promoter=self.sg_promoter,edit=self.ala_codon)
                oligos.append(oligo)
        else:
            slugs = []
            for i,tgt in enumerate(self.targets):
                slugs.extend(self.build_slug(tgt,n=3))

            for slug in slugs:
                if self.ha_margins=='max':
                    # spp + HA/2 + slug_len + HA/2 + sg_promoter + crRNA + end_primer = max_len
                    len_static = len(sp_primer) + len(slug) +len(self.sg_promoter) + crispr_rna_len + 20
                    ha_margins = (max_len - len_static)/2
                else:
                    ha_margins=self.ha_margins

                self.cassette_args = dict(
                    slug=slug,
                    sg_promoter=self.sg_promoter,
                    # id=self.lib_name+'.'+str(i),
                    name=self.lib_name,
                    description=self.lib_name,
                    up_margin=int(ha_margins),
                    down_margin=int(ha_margins),
                    crispr_len=crispr_rna_len
                )
                ec = self.build_cassette(**self.cassette_args)

                oligo = ec.assemble_oligo(sp_primer=sp_primer,edit=self.ala_codon)
                oligos.append(oligo)

        return oligos
    
# TODO: class CodonSubstitution()
# TODO: class TruncationScan()
# TODO: class SingleBaseDeletion()
# TODO: class DoubleBaseDeletion()
# TODO: class CodonDeletionScan()