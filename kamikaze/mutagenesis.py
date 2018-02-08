from .components import CassetteFactory

class AlanineScan(CassetteFactory):
    """
    Generates an alanine scan 

    ...

    Attributes
    ----------
    region  : slice
        Slice object representing CDS of sequence
    targets : list
        List of Slice objects for each codon

    Methods
    -------
    gen_payloads()
        Returns a list of payloads for each target in self.targets
    gen_edit_cassettes()
        Returns a list of edit cassettes for each target in self.targets
    
    build_slug(target,pam_site=None)
        Build slug for target and pam_site pair.
        Will find nearest pam_site to target if not defined

    """
    def __init__(self,*args):
        super().__init__(*args)
        targets = []
        i = self.region.start
        while i+3 < self.region.stop:
            targets.append(slice(i,i+3))
            i+=3
        self.targets = targets 
    
    def gen_payloads(self):
        return [self.build_payload(tgt) for tgt in self.targets]

    def gen_edit_cassettes(self):
        return [self.build_edit_cassette(tgt) for tgt in self.targets]

# TODO: class CodonSubstitution()
# TODO: class TruncationScan()
# TODO: class SingleBaseDeletion()
# TODO: class DoubleBaseDeletion()
# TODO: class CodonDeletionScan()