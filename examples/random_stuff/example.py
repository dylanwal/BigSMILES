import logging

import bigsmiles

# import chemdraw

# bigsmiles.Config.color_output = True

cases = [
    "F/C(C)=C(C(OC(NBr)C)C(C)(F)C)\C(OC(N)C)C(F)=C",

    "F/C=C/F",  # E-difluoroethene
    "F\C=C\F",  # Z-difluoroethene
    "C(\F)=C/F",  # trans
    "F\C=C/F",  # cis
    "F/C=C\F",  # cis
    "F/C(CC)=C/F",
    "F/C=C=C=C/F",  # trans
    "F/C=C=C=C\F", # cis
    "F/C=C/C/C=C\C",
    "F/C=C/CC=CC",  # partially specified
    "CC(F)/C=C/F",
    "CCC(\F)=C/F",  #
    "CCC(\F)=C(/F)CC",
    "CC/C(\F)=C/F",  # error
    "C(/[Ru])=C\C",
    "CC/C(F)=N\F",

    # "{[][<]OCCO[<],[>]C(=O)c1ccc(cc1)C(=O)[>],[>][H],[<]O[]}",
    # "C(={[>][<]=CC=[>][<]}=1)CCCCCC1",
    # "C=1CCC(={[>][<]=CC=[>][<]}=1)CCC",
    # "C({[>][<]CCO[>][<]}1)CCCCCC1",
    # "{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}",
    # "CCC1(CO)C(O)CC(O1)N",
    # "F/C=C/F",
    # "C1=CC=CC=C1C",
    # "F{[$][$]CC(CC[$])[$][$]}O",
    # "OC{[>][<]C(=O)OCC(=O)[<],[>]NCCCCCC(C)N[>][<]}CF",
    # "CC{[<]CC(C[>])C(C[<])CO[>]}O",
    # "CC(F){[<1][>1]NN[<1],[>1]CC(C)[<2],[>2]OO[<2][>2]}CCO",
    # "F{[<][>]CC(F)[<],[>]CCO[<][>]}P",
    # "CC{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}O",
    # "F{[<][>]OO[>],[>]C(N[<])C[<][>]}S",
]


def single():
    logging.basicConfig(level=logging.DEBUG)

    case = cases[0]
    print("case: ", case)
    output = bigsmiles.BigSMILES(case)
    print("output: ", output)
    for i, bond in enumerate(output.bonds):
        print(i, "   ", bond, bond.bond_order, bond.double_bond_ez)
    # output.print_tree(print_repr=False)

    # drawer = chemdraw.Drawer(str(output))
    # drawer.draw_html(auto_open=True)


def multiple():
    results = []
    for i, case in enumerate(cases):
        output = bigsmiles.BigSMILES(case)
        print(i, "  ", case, " --> ", output)
        results.append(str(output))
        for bond in output.bonds:
            print(bond, bond.bond_order, bond.double_bond_ez)
    # drawer = chemdraw.GridDrawer(results)
    # drawer.draw_html(auto_open=True)


if __name__ == "__main__":
    single()
    # multiple()


