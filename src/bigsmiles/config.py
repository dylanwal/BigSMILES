
class TerminalColors:
    Bold = '\033[1m'
    Underline = '\033[4m'
    Italics = '\33[3m'

    Red = '\33[91m'
    Green = '\33[92m'
    Yellow = '\33[93m'
    Blue = '\33[94m'
    Purple = '\33[95m'
    Cyan = '\33[96m'
    White = '\33[97m'
    Black = '\33[90m'

    Clear = '\33[0m'


class Config:
    color_output = False
    show_bond_descriptor_zero_index = False

    elements = {
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
        'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
        'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
        'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
        'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og', 'D'
    }

    elements_ordered = sorted(elements, reverse=True)  # sorted in reverse so "Cl" hits before "C"

    aromatic = {
        "b",
        "c",
        "n",
        "o",
        "p",
        "s"
    }

    elements_aromatic = elements.union(aromatic)

    organics = {
        "B",
        "C",
        "N",
        "O",
        "P",
        "S",
        "F",
        "Cl",
        "Br",
        "I"
    }

    atom_valences = {
        # put smallest valance first in tuple
        "H": (1,),
        "D": (1,),
        "Li": (1,),
        "Be": (2,),
        "B": (3,),
        "C": (4,),
        "N": (3,),
        "O": (2,),
        "F": (1,),
        "Na": (1,),
        "Mg": (2,),
        "Si": (4,),
        "P": (3, 5),
        "P5": (5,),  #
        "S": (2, 4, 6),
        "S4": (4,),  #
        "S6": (6,),  #
        "Cl": (1, 3, 5, 7),
        "K": (1,),
        "Ca": (2,),
        "Ga": (3,),
        "Ge": (4,),
        "As": (3, 5),
        "Se": (2, 4, 6),
        "Br": (1,),
        "Rb": (1,),
        "Sr": (2,),
        "In": (3,),
        "Sn": (2, 4),
        "Sb": (3, 5),
        "Te": (2, 4, 6),
        "I": (1, 3, 5, 7),
        "Cs": (1,),
        "Ba": (2,),
        "Ti": (3,),
        "Pb": (2, 4),
        "Bi": (3, 5),
        "Po": (2, 4, 6),
        "At": (1, 3, 5, 7),
        "Fr": (1,),
        "Ra": (2,),
    }

    @staticmethod
    def get_atom_possible_valence(symbol: str) -> tuple:
        try:
            return Config.atom_valences[symbol]
        except KeyError:
            return (8,)

    @classmethod
    def add_color(cls, text: str, color: str, skip_color: bool = False) -> str:
        if cls.color_output and not skip_color:
            return getattr(TerminalColors, color) + text + TerminalColors.Clear

        return text
