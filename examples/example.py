import bigsmiles

import chemdraw

# bigsmiles.Config.color_output = True

cases = [
    "N1=NN=C[N]1",  # None for Hydrogen
    "C(={[>][<]=CC=[>][<]}=1)CCCCCC1",
    "C(={[>][<]=CC=[>][<]}1)CCCCCC=1",
    "C(={[>][<]=CC=[>][<]}=1)CCCCCC=1",
    "C1CCC(={[>][<]=CC=[>][<]}=1)CCC",
    "C=1CCC(={[>][<]=CC=[>][<]}1)CCC",
    "C=1CCC(={[>][<]=CC=[>][<]}=1)CCC",

    "C({[>][<]CCO[>][<]}1)CCCCCC1",
    "C({[>][<]CCO[>][<]}12)CCCCCC12",

    "{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}",
    "CCC1(CO)C(O)CC(O1)N",
    "F/C=C/F",
    "C1=CC=CC=C1C",
    "F{[$][$]CC(CC[$])[$][$]}O",
    "OC{[>][<]C(=O)OCC(=O)[<],[>]NCCCCCC(C)N[>][<]}CF",
    "CC{[<]CC(C[>])C(C[<])CO[>]}O",
    "CC(F){[<1][>1]NN[<1],[>1]CC(C)[<2],[>2]OO[<2][>2]}CCO",
    "F{[<][>]CC(F)[<],[>]CCO[<][>]}P",
    "CC{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}O",
    "F{[<][>]OO[>],[>]C(N[<])C[<][>]}S",
]


def single():
    case = cases[0]
    print("case: ", case)
    output = bigsmiles.BigSMILES(case)
    print("output: ", output)
    output.print_tree(print_repr=False)

    # drawer = chemdraw.Drawer(str(output))
    # drawer.draw_html(auto_open=True)


def multiple():
    results = []
    for i, case in enumerate(cases):
        output = bigsmiles.BigSMILES(case)
        print(i, "  ", case, " --> ", output)
        results.append(str(output))

    # drawer = chemdraw.GridDrawer(results)
    # drawer.draw_html(auto_open=True)


if __name__ == "__main__":
    single()
    # multiple()
