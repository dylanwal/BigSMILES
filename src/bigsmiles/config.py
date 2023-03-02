
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

    @classmethod
    def add_color(cls, text: str, color: str, skip_color: bool = False) -> str:
        if cls.color_output and not skip_color:
            return getattr(TerminalColors, color) + text + TerminalColors.Clear

        return text
