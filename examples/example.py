
import bigsmiles

# bigsmiles.Config.color_output = True


def main():
    polymer_string = "CC{[][$]CC[]}"
    polymer = bigsmiles.BigSMILES(polymer_string)

    print("input: ", polymer_string)
    polymer.print_tree(print_repr=True)

    print("hi")


if __name__ == "__main__":
    main()
