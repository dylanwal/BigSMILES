
import bigsmiles

# bigsmiles.Config.color_output = True


def main():
    polymer_string = "{[][$]CC[$],[$]CC(CC)[$][]}"
    polymer = bigsmiles.BigSMILES(polymer_string)

    polymer.print_tree()


if __name__ == "__main__":
    main()
