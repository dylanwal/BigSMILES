import re

import bigsmiles.errors as errors
from bigsmiles.data_structures.bigsmiles import BigSMILES


rxn_split_pattern = re.compile("(?<!\[)(>>|>)")


def parse_reaction(text: str) -> tuple[list[BigSMILES], list[BigSMILES], list[BigSMILES]]:
    text = text.replace(" ", "")
    result = re.split(rxn_split_pattern, text)

    if ">>" in result:
        if ">" in result:
            raise errors.TokenizeError("Invalid BigSMILES reaction. Reactions must follow one of these patterns: "
                                       "\n\t 'reactants >> products'"
                                       "\n\t 'reactants > agents > products' ")

        # the pattern should be [reactants, >>, products]
        if result[1] != ">>":
            raise errors.TokenizeError("Invalid BigSMILES reaction. More than one '>>' detected.")

        reactants = process_chemical_block(result[0])
        products = process_chemical_block(result[2])

        return reactants, [], products

    # the pattern should be [reactants, >, agents, >, products]
    if result[1] != ">" and result[3] != ">":
        raise errors.TokenizeError("Invalid BigSMILES reaction. "
                                   "The pattern should be: 'reactants > agents > products'.")

    reactants = process_chemical_block(result[0])
    agents = process_chemical_block(result[2])
    products = process_chemical_block(result[4])

    return reactants, agents, products


def process_chemical_block(text: str) -> list[BigSMILES]:
    text_list = text.split(",")

    chemicals = []
    for chunk in text_list:
        chemical = BigSMILES(chunk)
        chemical = split_chemical(chemical)
        chemicals += chemical

    return chemicals


def split_chemical(bigsmile_: BigSMILES) -> list[BigSMILES]:
    if not bigsmile_.has_disconnect:
        return [bigsmile_]

    # check if

    return [bigsmile_]


def run_local():
    cases = [
        "C=CCBr >> C=CCI",
        "[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI",
        "C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]",
        "CC(=O)[OH]>>CC(=O)OCC",
        "CC(=O)[OH] >> CC(=O)OCC",  # extra spaces
        "C=Cc1ccccc1.C[CH2-]([Li+])CC"
        "C=C(C(=O)OC),CCCCSC(=S)SC(C)C(=O)O>CS(=O)C>{[>][<]CC(C(=O)OC)[>][<]}",  # RAFT
        "[CH2:1]=[CH:2][CH:3]=[CH:4][CH2:5][H:6]>> [H:6][CH2:1][CH:2]=[CH:3][CH:4]=[CH2:5]",
        "CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC",
    ]

    for case in cases:
        case = case.replace(" ", "")
        output = parse_reaction(case)
        print(output)


if __name__ == "__main__":
    run_local()
