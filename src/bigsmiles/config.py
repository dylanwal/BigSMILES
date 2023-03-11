
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
    # Will show colors in terminal when any object is printed
    color_output = False

    # True: [$1]; False: [$]
    show_bond_descriptor_one_index = False

    # True: will add ':' symbol for aromatic_elements bonds, False: hides the ':' in string outputs
    show_aromatic_bond = True

    # for small molecules only. True: 'C=1CCCCC=1'  False: 'C=1CCCCC1'
    show_multi_bonds_on_both_ring_index = False

    @classmethod
    def add_color(cls, text: str, color: str, skip_color: bool = False) -> str:
        if cls.color_output and not skip_color:
            return getattr(TerminalColors, color) + text + TerminalColors.Clear

        return text
