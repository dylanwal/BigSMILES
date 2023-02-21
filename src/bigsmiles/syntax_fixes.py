
from bigsmiles.bigsmiles import BigSMILES, Branch, StochasticObject, StochasticFragment


def run_syntax_fixes(bigsmiles: BigSMILES):
    if not bigsmiles:
        return

    remove_unnecessary_branch_symbols(bigsmiles)
    check_ring_numbers(bigsmiles)


def remove_unnecessary_branch_symbols(obj: BigSMILES | StochasticObject | StochasticFragment | Branch):
    # check if branch is notation is need; if not needed remove it and re-do ids
    # Example: CCC(C(CC)) --> CCCCCC

    for node in obj.nodes:
        if hasattr(node, "nodes"):
            remove_unnecessary_branch_symbols(node)

    branch_check(obj)


def branch_check(obj: StochasticObject | StochasticFragment | Branch):
    # checks if branch is end of nodes; if it is then it's not needed
    if isinstance(obj.nodes[-1], Branch):
        branch = obj.nodes.pop()
        obj.nodes += branch.nodes

        if hasattr(obj, 'parent'):
            renumber_branches(obj.parent)
        else:
            renumber_branches(obj)


def renumber_branches(obj, counter: int = 0) -> int:
    # Recursive function: loops through obj and re-numbers branch.id_
    for node in obj.nodes:
        if isinstance(node, Branch):
            node.id_ = counter
            counter += 1
        if hasattr(node, "nodes"):
            counter = renumber_branches(node, counter)

    return counter


def check_ring_numbers(obj: BigSMILES):
    # check if ring_counter had been skipped and restart everything at 1
    for i, ring in enumerate(obj.rings):
        ring.ring_id = i + 1
