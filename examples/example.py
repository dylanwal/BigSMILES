
import bigsmiles

bigsmiles.Config.color_output = True


def main():
    polymer_string = "CC(CC){[<1][>1]CC(C)[<2][>2]}CCO"  # "CC{[>][<]CC(C)[>][<]}CC(C)=C"
    polymer_string = "O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"

    polymer = bigsmiles.BigSMILES(polymer_string)
    polymer.nodes[4].to_string()
    print("input: ", polymer_string)
    polymer.print_tree(print_repr=True)

    g = polymer.graph()
    print(g)
    g.draw()


if __name__ == "__main__":
    main()
