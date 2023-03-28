
import bigsmiles


poly = bigsmiles.BigSMILES("[H]{[>][<]CC(C1=CC=CC=C1)[>][<]}[H]")

# BigSMILES
print(poly)
print(poly.atoms)
print(poly.bonds)
print(poly.nodes)
print(poly.molecular_formula)
print(poly.molar_mass)
print("")


# Atom
atom = poly.atoms[2]
print(atom)
print(atom.details)
print(atom.bonds)
print(atom.implicit_hydrogens)
print("")

# Bond
bond = poly.bonds[3]
print(bond)
print(bond.details)
print(bond.atom1, "-->", bond.atom2)
print("")

