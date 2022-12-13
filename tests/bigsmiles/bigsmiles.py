def local_test():
    # test_str = "CC1=C(C(=NC=C1C=O)C)O"
    test_str = "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"

    Config.color_output = True
    result = BigSMILES(test_str)
    print(test_str)
    print(result)
    print(repr(result))
    print()
    result.print_tree(print_repr=True)
    print("done")


if __name__ == "__main__":
    local_test()
