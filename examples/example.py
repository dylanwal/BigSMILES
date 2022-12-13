
import bigsmiles


def main():
    polymer_string = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
    polymer = bigsmiles.BigSMILES(polymer_string)
    polymer.print_tree()


if __name__ == "__main__":
    main()
