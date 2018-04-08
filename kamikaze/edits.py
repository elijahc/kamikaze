class CodonEdit():
    def __init__(self,codon):
        self.codon = codon

class Swap(CodonEdit):
    def __init__(self,codon,new_codon):
        super().__init__(codon)
        self.new_codon = new_codon

    def __call__(self):
        return self.new_codon

    def __str__(self):
        return str(self.codon) + ' to ' + str(self.new_codon)



"""
    def swap(self,new_codon='GCU'):
        return new_codon

    def trunc(self):
        return 'TAA'

    def delete(self):
        return ''
"""
