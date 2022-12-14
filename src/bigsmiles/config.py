
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

    atoms_with_valence = {
        # put smallest valance first in tuple
        "H": {"valence": (1,)},
        "D": {"valence": (1,)},
        "B": {"valence": (3,)},
        "C": {"valence": (4,)},
        "N": {"valence": (3,)},
        "O": {"valence": (2,)},
        "F": {"valence": (1,)},
        "Si": {"valence": (4,)},
        "P": {"valence": (3, 5)},
        "S": {"valence": (2, 4, 6)},
        "Cl": {"valence": (1, 3, 5, 7)},
        "Br": {"valence": (1,)},
        "I": {"valence": (1, 3, 5, 7)},
        "As": {"valence": (3, 5)},
    }

    @staticmethod
    def colors(text: str, color: str) -> str:
        return getattr(TerminalColors, color) + text + TerminalColors.Clear
