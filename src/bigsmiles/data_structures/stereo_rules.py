"""

Code for determining double bond E/Z stereo-chemistry

"""

from __future__ import annotations

import enum
import logging

import bigsmiles.reference_data.chemical_data as chemical_data
from bigsmiles.data_structures.bigsmiles import Bond, Atom, StochasticObject, BondDescriptorAtom


MAX_DEPTH = 100  # 100 is the max depth of traversal
DEFAULT_SCORE = 150


def get_other_node(bond: Bond, not_this_node: Atom | StochasticObject) \
        -> Atom | StochasticObject | BondDescriptorAtom | None:
    for atom in bond:
        if atom is not not_this_node:
            return atom

    return None


def get_other_bonds(atom: Atom, not_this_bond: tuple[Bond, ...] | list[Bond, ...]) -> list[Bond, ...]:
    bonds = []
    for bond in atom.bonds:
        if bond not in not_this_bond:
            bonds.append(bond)

    return bonds


def get_bond_with_ez_symbol(atom: Atom) -> Bond | None:
    ez_bond = None
    for bond in atom.bonds:
        if bond.symbol == "/" or bond.symbol == "\\":
            if ez_bond is not None:
                logging.warning(r"Extra '\' or '\\' detected." +
                                f"\nIgnoring symbol: {ez_bond.details} \nBigSMILES: {atom.root}")

            ez_bond = bond

    return ez_bond


class PriorityScore:
    __slots__ = ["atomic_number", "next_atom_highest_atomic_number", "number_bonds"]

    def __init__(self, atomic_number: int, next_atom_highest_atomic_number: int, number_bonds: float | int):
        self.atomic_number = atomic_number
        self.next_atom_highest_atomic_number = next_atom_highest_atomic_number
        self.number_bonds = number_bonds

    def __str__(self):
        return ",".join((str(getattr(self, name)) for name in self.__slots__))

    def __eq__(self, other):
        if self.atomic_number != other.atomic_number or self.number_bonds != other.number_bonds:
            return False

        return True

    def __lt__(self, other):
        if self.atomic_number < other.atomic_number:
            return True

        if self.atomic_number == other.atomic_number:
            if self.next_atom_highest_atomic_number < other.next_atom_highest_atomic_number:
                return True

            if self.next_atom_highest_atomic_number == other.next_atom_highest_atomic_number:
                if self.number_bonds < other.number_bonds:
                    return True

        return False


bond_scores = {
    ".": 0,
    "": 1,
    "-": 1,
    "/": 1,
    "\\": 1,
    ":": 1.5,
    "=": 2.1,  # little extra is to give them priority over two single bonds
    "#": 3.2,  # little extra is to give them priority over two single bonds
    "$": 4.1
}


def get_atom_score(node: Atom | StochasticObject | BondDescriptorAtom) -> int:
    if isinstance(node, Atom) and node.symbol in chemical_data.symbol_to_atomic_number:
        return chemical_data.symbol_to_atomic_number[node.symbol]

    return DEFAULT_SCORE


def get_score(node: Atom | StochasticObject | BondDescriptorAtom | None) -> PriorityScore:
    if node is None:
        return PriorityScore(0, 0, 0)

    # atom score
    atomic_number = get_atom_score(node)

    # next atom score
    next_atoms = (get_other_node(bond, node) for bond in node.bonds)
    next_atom_highest_atomic_number = max(get_atom_score(atom) for atom in next_atoms)

    # bond score
    if isinstance(node, Atom):
        number_bonds = sum(bond_scores[bond.symbol] for bond in node.bonds)
    elif isinstance(node, StochasticObject):
        number_bonds = node.bond_left.bond_order
    elif isinstance(node, BondDescriptorAtom):
        number_bonds = 0
    else:
        raise NotImplementedError("code bug")

    return PriorityScore(atomic_number, next_atom_highest_atomic_number, number_bonds)


class TraversalPath:

    def __init__(self, nodes: list = None, edges: list = None):
        self.nodes = nodes if nodes is not None else []
        self.edges = edges if edges is not None else []

    def __str__(self):
        return "->".join(repr(node) for node in self.nodes)

    def get_score(self) -> PriorityScore:
        return get_score(self.nodes[-1])


def next_step(path:  TraversalPath) -> list:
    additional_paths = []
    current_node = path.nodes[-1]

    if current_node is None:
        path.nodes.append(None)
        return additional_paths

    elif isinstance(current_node, Atom):
        atoms = []
        bonds = []
        for bond in current_node.bonds:
            if bond in path.edges:
                continue
            atoms.append(get_other_node(bond, current_node))
            bonds.append(bond)

        if len(atoms) == 0:
            path.nodes.append(None)
            return additional_paths

        atoms, bonds = sort_atoms_by_score(atoms, bonds)
        for i in range(len(atoms)):
            if i == 0:
                path.nodes.append(atoms[0])
                path.edges.append(bonds[0])
            else:
                additional_paths.append(
                    TraversalPath(path.nodes[:-1] + [atoms[i]], path.edges[:-1] + [bonds[i]])
                )

    elif isinstance(current_node, BondDescriptorAtom):
        path.nodes.append(None)

    elif isinstance(current_node, StochasticObject):
        if current_node.bond_left in path.edges:
            if current_node.bond_right is None:
                path.nodes.append(None)

            path.nodes.append(get_other_node(current_node.bond_right, current_node))
            path.edges.append(current_node.bond_right)

        else:
            if current_node.bond_left is None:
                path.nodes.append(None)

            path.nodes.append(get_other_node(current_node.bond_left, current_node))
            path.edges.append(current_node.bond_left)

    else:
        raise NotImplementedError("coding bug")

    return additional_paths


def sort_atoms_by_score(atoms: list[Atom, ...], bonds: list[Bond, ...]) -> tuple[list[Atom, ...], list[Bond, ...]]:
    score = []
    for atom in atoms:
        score.append(get_atom_score(atom))

    return [x for _, x in sorted(zip(score, atoms))], [x for _, x in sorted(zip(score, bonds))],


class TraversalDirection:

    def __init__(self, start_atom: Atom, bond_prior: Bond):
        self.start_atom = start_atom
        self.bond_prior = bond_prior
        self.paths: [TraversalPath] = [TraversalPath([start_atom], [bond_prior])]
        self.best_path = self.paths[0]
        self._scores = None
        self._up_to_date = False

    @property
    def num_paths(self) -> int:
        return len(self.paths)

    def scores(self) -> list[PriorityScore]:
        if not self._up_to_date:
            self._scores = [path.get_score() for path in self.paths]
            index_not_none = next(i for i, score in enumerate(self._scores) if score is not None)
            self.best_path = self.paths[index_not_none]
            self._up_to_date = True

        return self._scores

    def one_traverse_step(self):
        self._up_to_date = False

        # traverse the graph one step
        for i in range(len(self.paths)):  # index used as size of self.up may change during loop
            additional_paths = next_step(self.paths[i])
            self.paths += additional_paths


class TraversalOutcomes(enum.Enum):
    up = 1
    down = 2
    tie = 3


class TraversalSide:
    def __init__(self, double_bond: Bond, side: bool):
        """
        Parameters
        ----------
        double_bond:

        side:
            True: left, False: right
        """
        self.double_bond = double_bond
        self.side = side

        self.up: TraversalDirection | None = None
        self.down: TraversalDirection | None = None
        self.winner: TraversalOutcomes | None = None
        self.number_of_double_bonds_in_series: int = 0  # in a row after the starting double bond
        self.setup_successful: bool = True

        self.setup()

    def setup(self):
        if self.side:
            atom = self.double_bond.atom1
        else:
            atom = self.double_bond.atom2

        # get starting atom
        starting_atom = self.check_for_multiple_double_bonds(atom, self.double_bond)
        if starting_atom is None:
            return

        # determine up down directions
        if self.side:
            self.determine_up_down_atom_left(starting_atom)
        else:
            self.determine_up_down_atom_right(starting_atom)

    def check_for_multiple_double_bonds(self, atom: Atom, bond: Bond) -> Atom | None:
        """
        check for multiple double bonds in a row
        Maybe recursive

        Parameters
        ----------
        atom:
            current atom position
        bond:
            prior bond

        Returns
        -------
        atom:
            Atom to start up/down traversals from

        """
        other_bonds = get_other_bonds(atom, [bond])
        if not other_bonds:
            self.setup_successful = False
            return None  # double bond is terminated with implicit hydrogens or lone pairs
        if len(other_bonds) == 1 and other_bonds[0].symbol == "=":
            next_atom = get_other_node(other_bonds[0], atom)
            self.number_of_double_bonds_in_series += 1
            return self.check_for_multiple_double_bonds(next_atom, other_bonds[0])

        return atom

    def determine_up_down_atom_left(self, start_atom: Atom):
        bond_ez = get_bond_with_ez_symbol(start_atom)
        if bond_ez is None:
            self.setup_successful = False
            return  # no stereo bonds '/' or '\\'

        atom_ez = get_other_node(bond_ez, start_atom)

        # set up path
        direction = TraversalDirection(atom_ez, bond_ez)
        if bond_ez.symbol == "/":
            if atom_ez.parent is start_atom.parent:
                self.up = direction  # main chain
            else:
                self.down = direction  # branch
        else:  # "\\"
            if atom_ez.parent is start_atom.parent:
                self.down = direction  # main chain
            else:
                self.up = direction  # branch

        self.figure_out_other_bond(start_atom, bond_ez)

    def determine_up_down_atom_right(self, start_atom: Atom):
        bond_ez = get_bond_with_ez_symbol(start_atom)
        if bond_ez is None:
            self.setup_successful = False
            return  # no stereo bonds '/' or '\\'

        atom_ez = get_other_node(bond_ez, start_atom)

        # set first atom that has stereo symbol
        if bond_ez.symbol == "/":
            self.down = TraversalDirection(atom_ez, bond_ez)
        else:  # "\\"
            self.up = TraversalDirection(atom_ez, bond_ez)

        self.figure_out_other_bond(start_atom, bond_ez)

    def figure_out_other_bond(self, start_atom: Atom, bond_ez: Bond):
        # figure out other bond
        other_bond = get_other_bonds(start_atom, (self.double_bond, bond_ez))
        if len(other_bond) == 1:
            other_atom = get_other_node(other_bond[0], start_atom)
            if self.down is None:
                self.down = TraversalDirection(other_atom, other_bond[0])
            else:
                self.up = TraversalDirection(other_atom, other_bond[0])

        elif len(other_bond) == 0:
            # other bond must be implicit hydrogens or terminal
            if self.down is None:
                self.winner = TraversalOutcomes.up
            else:
                self.winner = TraversalOutcomes.down

        else:
            raise NotImplementedError("more than 2 bonds coming off double bond atom")

    def check_for_winner(self):
        if self.winner:
            return

        scores_up = self.up.scores()
        scores_down = self.down.scores()
        for i in range(max([len(scores_up), len(scores_down)])):
            if i > len(scores_up) - 1 and i > len(scores_up) - 1:
                self.winner = TraversalOutcomes.tie
            if i > len(scores_up) - 1:
                self.winner = TraversalOutcomes.down
                break
            if i > len(scores_down) - 1:
                self.winner = TraversalOutcomes.up
                break

            if scores_up[i] == scores_down[i]:
                if scores_up is None:
                    continue
                else:
                    break

            if scores_up > scores_down:
                self.winner = TraversalOutcomes.up
            else:
                self.winner = TraversalOutcomes.down

    def run_graph_traversal(self):
        self.check_for_winner()
        if self.winner:
            return

        logging.debug(f" Starting '{'left' if self.side else 'right'}' E/Z graph traversal:"
                      f" {self.double_bond.root.to_string(show_repr=(Atom,))}"
                      f"\n-------------------------------")
        logging.debug(f"\t {'step':>4}, {'up':<30}, {'num_up_paths':>12}, {'down':>30}, {'num_down_paths':>12}")
        logging.debug(f"\t {0:4}, {str(self.up.best_path):<30}, {self.up.num_paths:>12}"
                      f" {str(self.down.best_path):<30}, {self.down.num_paths:>12}")

        # perform graph traversal
        for i in range(1, MAX_DEPTH):
            self.up.one_traverse_step()
            self.down.one_traverse_step()

            self.check_for_winner()
            logging.debug(f"\t {i:4}, {str(self.up.best_path):<30}, {self.up.num_paths:>12}"
                          f" {str(self.down.best_path):<30}, {self.down.num_paths:>12}")

            if self.winner:
                return

        else:
            logging.warning(f'Max traversal depth reached on E/Z determination of: '
                            f'\nBigSMILES: {self.double_bond.parent} '
                            f'\nBond: {repr(self.double_bond)}')
            return


def check_setup(left_side: TraversalSide, right_side: TraversalSide) -> bool:
    """ Check to see if there were reasons to terminate early from traversal setup """
    if not left_side.setup_successful or not right_side.setup_successful:
        return True

    # check for 'odd' number of double bonds in a row
    if (left_side.number_of_double_bonds_in_series + right_side.number_of_double_bonds_in_series + 1) % 2 == 0:
        # left + right + 1 (starting double bond) = even  --> return None for E/Z
        logging.warning("Double bond stereo-chemistry notation '/' and '\\' should "
                        f"not be used on even number of double bonds; they are axial chiral and use @/@@ notation on "
                        f"center atom.\nBigSMILES:{left_side.double_bond.parent}")
        return True

    return False


def determine_winner(left_side: TraversalSide, right_side: TraversalSide) -> str | None:
    if left_side.winner is None or right_side.winner is None:
        return None

    if left_side.winner is TraversalOutcomes.tie and right_side.winner is TraversalOutcomes.tie:
        logging.warning("Double bond stereo-chemistry notation '/' and '\\' should "
                        "not be used on symmetric molecules.")
        return None

    if (left_side.winner is TraversalOutcomes.tie or right_side.winner is TraversalOutcomes.tie) or \
            (left_side.winner == right_side.winner):
        return "Z"

    return "E"


def get_double_bond_ez(double_bond: Bond) -> str | None:
    """
    Determines the E/Z stereochemistry of a double Bond

    Parameters
    ----------
    double_bond:
        double bond to start analysis at

    Returns
    -------
    result:
        will return 'E'/'Z' for odd number of double bonds, None for even number of double bonds (axial chiral)

    Notes
    -----

    * Breath-first-search (BFS) prioritizing heavy atoms and following more bonds.
    Starts out at the double bond with two traverses heading left and two right,
    with one going up and the other going down for a total of 4 simultaneous BFS.
    At each step the below rules are applied to compute a score
    The searches terminate once a difference between the up-down score.
    * Rules:
        * rule 1: higher atomic mass has higher priority
        * rule 2: (if same atomic mass) atom with more of these substituent
        * rule 3: propagate down chain till there is a difference
        * double/triple bonds are considered to be bonded to an equivalent single bonds C=C --> C(C)-C(C)
    * stochastic objects are scored as 150 and treated as a single atom
    * max traversal depth is 200

    """
    left_side = TraversalSide(double_bond, True)
    right_side = TraversalSide(double_bond, False)

    if check_setup(left_side, right_side):
        return None

    left_side.run_graph_traversal()
    right_side.run_graph_traversal()

    return determine_winner(left_side, right_side)
