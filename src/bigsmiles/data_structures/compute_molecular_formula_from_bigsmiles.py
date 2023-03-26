
from bigsmiles.data_structures.bigsmiles import BigSMILES, StochasticFragment, Atom, StochasticObject, Branch
from bigsmiles.data_structures.molecular_formula import MolecularFormula


def compute_molecular_formula_from_bigsmiles(bigsmiles: BigSMILES | StochasticFragment) -> MolecularFormula:
    element_dict = _compute_molecular_formula_from_bigsmiles_loop(bigsmiles)
    return MolecularFormula(element_dict)


def _compute_molecular_formula_from_bigsmiles_loop(obj: BigSMILES | StochasticFragment | Branch) -> dict[str, int]:
    element_dict: dict[str, int] = {"H": 0}

    for node in obj.nodes:
        if isinstance(node, Atom):
            add_atom_to_dict(element_dict, node)
        elif isinstance(node, StochasticObject):
            add_stochastic_object_to_dict(element_dict, node)
        elif isinstance(node, Branch):
            merge_dict(element_dict, _compute_molecular_formula_from_bigsmiles_loop(node))
        # skip Bond, Bonding descriptor atom; stochastic fragment process separately

    if element_dict["H"] == 0:
        del element_dict["H"]

    return element_dict


def add_atom_to_dict(element_dict: dict[str, int], atom: Atom):
    if atom.symbol in element_dict:
        element_dict[atom.symbol] += 1
    else:
        element_dict[atom.symbol] = 1

    element_dict["H"] += int(atom.implicit_hydrogens)


def add_stochastic_object_to_dict(element_dict: dict[str, int], stochastic_object: StochasticObject):
    if "{}" in element_dict:
        element_dict["{}"] += 1
    else:
        element_dict["{}"] = 1

    # recursive compute all molecular formulas
    for stochastic_fragment in stochastic_object.nodes:
        stochastic_fragment._molecular_formula = compute_molecular_formula_from_bigsmiles(stochastic_fragment)


def merge_dict(dict1: dict[str, int], dict2: dict[str, int]):
    """ merge dict2 into dict1 """
    for k, v in dict2.items():
        if k in dict1:
            dict1[k] += dict2[k]
        else:
            dict1[k] = dict2[k]
