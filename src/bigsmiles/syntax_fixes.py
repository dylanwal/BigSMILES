
from bigsmiles.bigsmiles import BigSMILES, Branch, Atom, Bond

def remove_unnecessary_branch_symbols():
    # check if branch is notation is need; if so remove it
    # Example: CCC(C(CC)) --> CCC(CCC)


def ring_notation(bigsmiles: bigsmiles):
    ## C12 -> C=1 or C123 -> C#1

    pairs = set()
    for ring in bigsmiles.rings:
        if

