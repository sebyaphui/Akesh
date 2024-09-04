"""
from simulation import Phenotype, SequencePhenotypeStructure

struct = SequencePhenotypeStructure([])
p = Phenotype(struct, 2)
p2 = Phenotype(struct, 2)

a = {}


a[p] = 2
a[p2] = 4

print(a[p])
print(a[p2])
"""

from discrete_model import Phenotype, LatticePhenotypeStructure, SequencePhenotypeStructure


"""
phen_struct = LatticePhenotypeStructure(1, 100)
phen_1 = Phenotype(phen_struct, 10)

a = {}

a[phen_1] = 1


phen_2 = Phenotype(phen_struct, 10)

print(phen_2 in a)
"""



##peptides = read_phenotypes("peptides_from_TULIP.txt")
#print(peptides)