
import bigsmiles

bigsmiles.Config.color_output = True


def main():
    polymer_string = "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"
    polymer = bigsmiles.BigSMILES(polymer_string)

    print("input: ", polymer_string)
    polymer.print_tree(print_repr=False)


if __name__ == "__main__":
    main()
