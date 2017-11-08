from .components import CassetteFactory

class AlanineScan(CassetteFactory):
    
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

    def gen_edit_cassettes(self,mut='GCU'):
        return [self.build_edit_cassette(tgt) for tgt in self.targets]