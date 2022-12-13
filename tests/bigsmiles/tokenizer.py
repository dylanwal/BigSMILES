def local_test():
    """ For quick testing. """
    # test_str = "C#[12C](C12O)Ncsn[C@][C@@][Fe+][Fe++][Fe+2][CH2]Ru%12CCC/C=C/C[OH3+]"
    # test_str = "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"
    test_str = "{[]CCC([$1[$1]1])CC(OC[$1[$2]1])(CC[$1[$1]2])OC[$1[$2]2][]}"
    result = tokenize(test_str)
    for r in result:
        print(r)