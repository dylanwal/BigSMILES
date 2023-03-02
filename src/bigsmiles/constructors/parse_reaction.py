import re

rxn_split_pattern = "(>>|>)"


def parse_reaction(text: str):
    text = text.strip()

    result = re.split(rxn_split_pattern, text)

    return result


def run_local():
    cases = [
        "C=CCBr>>C=CCI",
        "[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI",
        "C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]",
        "CC(=O)[OH]>>CC(=O)OCC",
        "CC(=O)[OH] >> CC(=O)OCC",  # extra spaces
        "C=Cc1ccccc1.C[CH2-]([Li+])CC"
        "C=C(C(=O)OC),CCCCSC(=S)SC(C)C(=O)O>CS(=O)C>{[>][<]CC(C(=O)OC)[>][<]}",  # RAFT
        "[CH2:1]=[CH:2][CH:3]=[CH:4][CH2:5][H:6]>> [H:6][CH2:1][CH:2]=[CH:3][CH:4]=[CH2:5]",
        "CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC"
    ]

    for case in cases:
        output = parse_reaction(case)
        print(output)


if __name__ == "__main__":
    run_local()
