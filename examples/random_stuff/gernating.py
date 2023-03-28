
import bigsmiles


# use defaults
polymer = bigsmiles.Polymer("[H]{[>][<]CC(C1=CC=CC=C1)[>][<]}[H]")
molecule = bigsmiles.methods.generate_molecules(polymer)
print(molecule)


# define manually
polymer = bigsmiles.Polymer("[H]{[>][<]CC(C1=CC=CC=C1)[>][<]}[H]")
polymer.spec[0].distribution = bigsmiles.distributions.LogNormal(Mn=12_000, D=1.12)
molecule = bigsmiles.methods.generate_molecules(polymer)
print(molecule)


# define with string
polymer = bigsmiles.Polymer("{[][$]C([$])C=O,[$]CC([$])CO;[$][H], [$]O[]}|flory_schulz(0.0011)|")
molecule = bigsmiles.methods.generate_molecules(polymer)
print(molecule)

"""

linear 
* homopolymer
    * N
* diblock (two stochastic objects)
    * N for each block
* diblock (bonding descriptor index)
    * N for each block (index)
* random
    * N
    * % of each 
* gradient 
    * N 
    * function that return % chance of each monomer based on length down backbone
* alternating 
    * N
* rings
    * N   
    
graft
* bottlebrush
    * N backbone
    * N brush;  func(N_backbone, N_brush)
    
star
* N of each arm

multiple endgroups
* percent on end groups

hyperbranch/dendramers
* N
* % each (function that return % chance of each monomer based on length down each branch)

* implicit end groups terminate with [H] or *


no ladder; no networks
"""