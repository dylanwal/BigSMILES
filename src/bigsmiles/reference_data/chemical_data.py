import os
import json


bond_mapping = {
    None: None,
    ".": 0,
    "": 1,
    "-": 1,
    "/": 1,
    "\\": 1,
    ":": 1.5,
    "=": 2,
    "#": 3,
    "$": 4
}

stereo_bonds = {"/", "\\"}

# load element data from JSON
file_path = os.path.dirname(os.path.realpath(__file__))
path_to_data = os.path.join(file_path, "elements_data.json")
with open(path_to_data, "r", encoding="UTF-8") as f:
    element_data = json.load(f)


class Element:
    __slots__ = ("name", "atomic_number", "atomic_mass", "period", "symbol", "electron_configuration", "group",
                 "aromatic", "organics", "conflict", "valences")

    def __init__(self,
                 name: str,
                 atomic_number: int,
                 atomic_mass: float,
                 period: int,
                 symbol: str,
                 electron_configuration: str,
                 group: str,
                 aromatic: bool,
                 organics: bool,
                 conflict: bool,
                 valences: tuple[int]
                 ):
        self.name = name
        self.atomic_number = atomic_number
        self.atomic_mass = atomic_mass
        self.period = period
        self.symbol = symbol
        self.electron_configuration = electron_configuration
        self.group = group
        self.aromatic = aromatic
        self.organics = organics
        self.conflict = conflict
        self.valences = valences

    def __str__(self):
        return self.symbol

    def __repr__(self):
        return self.symbol


periodic_table = tuple(Element(atomic_number=k, **v) for k, v in element_data.items())
del element_data

# create tuples, sets, and dictionaries for performance
element_symbols = {element.symbol for element in periodic_table}
elements_ordered = sorted(element_symbols, reverse=True)  # sorted in reverse so "Cl" hits before "C"
aromatic_elements = {element.symbol.lower() for element in periodic_table if element.aromatic}
elements_aromatic = element_symbols.union(aromatic_elements)
organic_elements = {element.symbol for element in periodic_table if element.organics}
organic_ordered = sorted(organic_elements, reverse=True)  # sorted in reverse so "Cl" hits before "C"
atom_valences = {element.symbol: element.valences for element in periodic_table}
