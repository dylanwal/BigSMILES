
import bigsmiles

# bigsmiles.Config.color_output = True


def main():
    cases = [
    "C1=CC=CC=C1C",
    "F{[$][$]CC(CC[$])[$][$]}O",
    "OC{[>][<]C(=O)OCC(=O)[<],[>]NCCCCCC(C)N[>][<]}CF",
    "CC{[<]CC(C[>])C(C[<])CO[>]}O",
    "CC(F){[<1][>1]NN[<1],[>1]CC(C)[<2],[>2]OO[<2][>2]}CCO",
    "F{[<][>]CC(F)[<],[>]CCO[<][>]}P",
    "CC{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}O",
    "F{[<][>]OO[>],[>]C(N[<])C[<][>]}S",
    ]

    case = cases[0]
    # for case in cases:
    polymer = bigsmiles.BigSMILES(case)
    print("input: ", case)
    polymer.print_tree(print_repr=False)

    g = polymer.graph()
    print(g)
    g.draw_plotly()


if __name__ == "__main__":
    main()
