from __future__ import annotations
from bigsmiles.data_structures.bigsmiles import BigSMILES


class Reaction:
    """

    data structure for reaction BigSMILES

    !!! info "Accepted patterns"

        * reactant '>' agent '>' product
        * reactant '>>' product
        * reactant.reactant>>product.product
        * reactant,reactant>>product,product (comma supported for easier parsing)

    !!! info "Definitions"
        **^^reactants:^^** materials that contributing one or more atoms to the product

        **^^agents:^^** materials that don't contribute any atoms to the product or receive atoms from the reactant
        (catalysts, solvents)

        **^^product:^^** output of reactions. all atoms should have come from the reactants

    """
    __slots__ = ["reactants", "agents", "products", "__dict__"]

    def __init__(self, text: str | None = None, **kwargs):
        """

        Parameters
        ----------
        text: str
            Reaction BigSMILES string
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
        self.reactants: list[BigSMILES] = []
        self.agents: list[BigSMILES] = []
        self.products: list[BigSMILES] = []

        # parse input string
        if text:
            from bigsmiles.constructors.constructor_reaction import parse_reaction
            self.reactants, self.agents, self.products = parse_reaction(text)

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
