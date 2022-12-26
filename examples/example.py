
import bigsmiles

# bigsmiles.Config.color_output = True


def main():
    # polymer_string = "F{[$][$]CC(CC[$])[$][$]}O"
    # polymer_string = "OC{[>][<]C(=O)OCC(=O)[<],[>]NCCCCCC(C)N[>][<]}CF"
    # polymer_string = "CC{[<]CC(C[>])C(C[<])CO[>]}O"
    polymer_string = "CC(CC){[<1][>1]CC(C)[<2][>2]}CCO"

    polymer = bigsmiles.BigSMILES(polymer_string)

    print("input: ", polymer_string)
    polymer.print_tree(print_repr=True)

    g = polymer.graph()
    print(g)
    g.draw()


if __name__ == "__main__":
    main()
