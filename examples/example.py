
import bigsmiles

bigsmiles.Config.color_output = True


def main():
    polymer_string = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
    polymer = bigsmiles.BigSMILES(polymer_string)

    print("input: ", polymer_string)
    polymer.print_tree(print_repr=True)

    g = polymer.graph()
    print(g)
    g.draw()


if __name__ == "__main__":
    main()
