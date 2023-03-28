
import bigsmiles

polymer = bigsmiles.Polymer("[H]{[>][<]CC(C1=CC=CC=C1)[>][<]}[H]")

print(polymer)
# polymer.spec[0].distribution = bigsmiles.distribution.LogNormal(Mn=10_000, D=1.12)
#

molecule = bigsmiles.methods.generate_molecules(polymer)
print(molecule)

