from components import CassetteFactory

cf = CassetteFactory(fi='kras_mrna_va.fa',region=slice(192,761))

# Generate pam sites
# pam_sites = cf.pam_sites()

# Find nearest pam sites
shift = lambda slc,x: slice(slc.start+x,slc.stop+x)
tgt = shift(slice(192,195),63)

rev_sites = [ps for ps in cf.find_pam_sites() if ps.strand==-1]
rev_seqs = [ps.extract(cf.reference.seq) for ps in rev_sites]
pl = cf.build_payload(shift(tgt,9))
import ipdb; ipdb.set_trace()
