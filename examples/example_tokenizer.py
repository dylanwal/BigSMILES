
import bigsmiles

bigsmiles_str = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
bigsmiles_str = "{[][>]NCCCCCCN[>],[<]C(=O)CCCCC(=O)[<],[>]Cl,[<][H][]}"


# Text tokens
bigsmiles_split_str = bigsmiles.tokenize_text(bigsmiles_str)
print(bigsmiles_split_str)
# ['C', 'C', '{', '[>]', '[<]', 'C', 'C', '(', 'C', ')', '[>]', '[<]', '}', 'C', 'C', '(', 'C', ')', '=', 'C']

# Tokens
bigsmiles_tokens = bigsmiles.tokenize(bigsmiles_str)
print(bigsmiles_tokens)
# [
# Atom: C, Atom: C, StochasticStart: {, BondDescriptor: [>], BondDescriptor: [<], Atom: C, Atom: C,
# BranchStart: (, Atom: C, BranchEnd: ), BondDescriptor: [>], BondDescriptor: [<], StochasticEnd: },
# Atom: C, Atom: C, BranchStart: (, Atom: C, BranchEnd: ), Bond: =, Atom: C
# ]


# Atom
atom_dict = bigsmiles.tokenize_atom_symbol("[13C@H+:1]")
print(atom_dict)
# {'isotope': 13, 'symbol': 'C', 'stereo': '@', 'hydrogens': 1, 'charge': 1, 'class_': 1}


# Bonding Descriptor
bond_descr = bigsmiles.tokenize_bonding_descriptor("[$1]")
print(bond_descr)
# ('$', 1)
