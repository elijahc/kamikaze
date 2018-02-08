from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class PAM():
    """Base class for PAMs

    Attributes
    ----------
    strain : str
    seq : Seq
    """
    def __init__(self,strain=None, seq=None):
        if strain is not None and seq is None:
            self.strain = strain
            if strain.lower()=='cas9':
                self.seq = Seq('NGG',IUPAC.unambiguous_dna)

        elif strain is None and seq is not None:
            self.seq = seq
            if self.seq.upper()=='NGG': 
                self.strain='cas9'

class Cas9(PAM):
    """Cas9 PAM

    Attributes
    ----------
    strain : str
    seq : Seq
    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)