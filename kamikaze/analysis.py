import requests
import json

def wge_xref(guide,species):
    wge_endpoint = "http://www.sanger.ac.uk/htgt/wge/api/search_by_seq"
    if species not in ['Grch38','Grcm38']:
        raise "Species must be Grch38 or Grcm38"
    
    params = {'get_db_data':1,'pam_right':2,'seq':guide,'species':species}
    r = requests.get(wge_endpoint,params=params)
    print(r)
    return r.json()

class GuideAnalysis():
    def __init__(self,guides):
        guides = set(guides)
        self.guides = guides
    
    def wge(self,species):
        wge_results = []
        for g in self.guides:
            print('requesting guide '+str(g)+'...')
            r = wge_xref(g,species)
            wge_results.append(r)
        
        self.wge_results = wge_results

        return wge_results