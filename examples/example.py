import bigsmiles

import chemdraw

# bigsmiles.Config.color_output = True

cases = [
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
    for case in cases:
        output = bigsmiles.BigSMILES(case)
        print(case, " --> ", output)

    drawer = chemdraw.GridDrawer(results)
    drawer.draw_html(auto_open=True)


if __name__ == "__main__":
    single()
    # multiple()
