
from bigsmiles.data_structures.bigsmiles import BigSMILES


def extract_text_between_pipe(string: str) -> tuple[list[str, ...], list[str, ...]]:
    removed_text = []
    remaining_text = []
    while "|" in string:
        start = string.find("|")
        end = string.find("|", start + 1)
        removed_text.append(string[start+1:end])
        remaining_text.append(string[:start] + string[end+1:])
        string = string[:start] + string[end+1:]
    remaining_text.append(string)

    return removed_text, remaining_text


class Polymer:

    def __init__(self, text: str = None):
        if '|' in text:
            bigsmiles_text, gen_text = extract_text_between_pipe(text)
            text = "".join(bigsmiles_text)

        self.bigsmiles = BigSMILES(text)
        self.dis = None

    def get_BigSMILES_gen_str(self) -> str:
        pass

