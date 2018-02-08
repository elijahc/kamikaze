import os

package_dir = os.path.dirname(os.path.abspath(__file__))

def kras():
    return os.path.join(package_dir,'data','kras_mrna_va.fa')

def galk():
    return os.path.join(package_dir,'data','galk.fa')