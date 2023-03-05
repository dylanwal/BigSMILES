"""

Data structure for BigSMILE reaction notation.

reactant '>' agent '>' product
reactant '>>' product

reactant.reactant>>product.product
reactant,reactant>>product,product (comma supported for easier parsing)

agents: are not consumed during the reaction

"""

from bigsmiles.data_structures.bigsmiles import BigSMILES


class Reaction:
    __slots__ = ["reactants", "agents", "products", "__dict__"]

    def __init__(self, input_text: str = None, **kwargs):
        self.reactants: list[BigSMILES] = []
        self.agents: list[BigSMILES] = []
        self.products: list[BigSMILES] = []

        # parse input string
        if input_text:
            from bigsmiles.constructors.parse_reaction import parse_reaction
            self.reactants, self.agents, self.products = parse_reaction(input_text)

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        text = ""
        text += ",".join(chem.to_string(show_hydrogens, show_atom_index, print_repr, skip_color)
                         for chem in self.reactants)
        if self.agents:
            text += ">"
            text += ",".join(chem.to_string(show_hydrogens, show_atom_index, print_repr, skip_color)
                             for chem in self.agents)
            text += ">"
        else:
            text += ">>"
        text += ",".join(chem.to_string(show_hydrogens, show_atom_index, print_repr, skip_color)
                         for chem in self.products)

        return text
