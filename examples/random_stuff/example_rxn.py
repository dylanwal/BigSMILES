import bigsmiles

import chemdraw


cases = [
    "C=Cc1ccccc1.C[CH-](.[Li+])CC>Cc1ccccc1>CC(CC){[>][<]CC(c1ccccc1)[>][<]}[H]",

    "C=CCBr >> C=CCI",
    "[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI",
    "C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]",
    "CC(=O)[OH]>>CC(=O)OCC",
    "CC(=O)[OH] >> CC(=O)OCC",  # extra spaces
    "C=Cc1ccccc1.C[CH-](.[Li+])CC>Cc1ccccc1>CC(CC){[>][<]CC(c1ccccc1)[>][<]}[H]",
    "C=C(C(=O)OC)>CCCCSC(=S)SC(C)C(=O)O,CS(=O)C>{[>][<]CC(C(=O)OC)[>][<]}",  # RAFT
    "[CH2:1]=[CH:2][CH:3]=[CH:4][CH2:5][H:6]>> [H:6][CH2:1][CH:2]=[CH:3][CH:4]=[CH2:5]",
    "CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC",
]


def single():
    case = cases[0]
    print("case: ", case)
    output = bigsmiles.Reaction(case)
    print("output: ", output)

    # drawer = chemdraw.Drawer(str(output))
    # drawer.draw_html(auto_open=True)


def multiple():
    results = []
    for i, case in enumerate(cases):
        output = bigsmiles.Reaction(case)
        print(i, "  ", case, " --> ", output)
        results.append(str(output))

    # drawer = chemdraw.GridDrawer(results)
    # drawer.draw_html(auto_open=True)


if __name__ == "__main__":
    single()
    # multiple()
