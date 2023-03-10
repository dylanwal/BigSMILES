"""

Code for determining double bond E/Z stereo-chemistry

"""

from __future__ import annotations

import enum
import logging

import bigsmiles.reference_data.chemical_data as chemical_data
from bigsmiles.data_structures.bigsmiles import Bond, Atom, StochasticObject, BondDescriptorAtom

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


def get_score(node: Atom | StochasticObject | BondDescriptorAtom | None) -> int:
    if node is None:
        return 0
    elif isinstance(node, Atom):
        return chemical_data.symbol_to_atomic_number[node.symbol]

    return DEFAULT_SCORE


class TraversalPath:
    def __init__(self,
                 current_node: Atom | StochasticObject | BondDescriptorAtom | None,
                 prior_edge: Bond | None
                 ):
        self.current_node = current_node
        self.prior_edge = prior_edge

    def next_step(self):
        if isinstance(self.current_node, BondDescriptorAtom):
            self.current_node = None

        elif isinstance(self.current_node, Atom):
            best_atom = None
            best_bond = None
            best_score = 0
            for bond in self.current_node.bonds:
                if bond is self.prior_edge:
                    continue
                atom = get_other_node(bond, self.current_node)
                score = get_score(atom)
                if score > best_score:
                    best_score = score
                    best_atom = atom
                    best_bond = bond

            self.prior_edge = best_bond
            self.current_node = best_atom

        elif isinstance(self.current_node, StochasticObject):
            if self.current_node.bond_left is self.prior_edge:
                bond = self.current_node.bond_right
                if bond is None:
                    self.prior_edge = None
                    self.current_node = None

                self.current_node = get_other_node(bond, self.current_node)
                self.prior_edge = bond
                return
            else:
                bond = self.current_node.bond_left
                if bond is None:
                    self.prior_edge = None
                    self.current_node = None

                self.current_node = get_other_node(bond, self.current_node)
                self.prior_edge = bond

        else:
            raise ValueError("coding bug")


class TraversalOutcomes(enum.Enum):
    up = 1
    down = 2
    tie = 3


class TraversalSide:
    def __init__(self, double_bond: Bond, side: bool):
        """

        Parameters
        ----------
        double_bond
        side:
            True: left, False: righ
        """
        self.double_bond = double_bond
        self.side = side

        self.up: TraversalPath | None = None
        self.down: TraversalPath | None = None
        self.winner: TraversalOutcomes | None = None
        self.number_of_double_bonds_in_series: int = 0  # in a row after the starting double bond
        self._seen_bond: set[int] = {double_bond.id_}
        self._seen_atom: set[int] = set()
        self.setup_successful: bool = True

        self.setup()

    @property
    def done(self) -> bool:
        if self.winner is None:
            return False
        return True

    def setup(self):
        if self.side:
            atom = self.double_bond.atom1
        else:
            atom = self.double_bond.atom2

        # get starting atom
        self._seen_atom.add(atom.id_)
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
            self._seen_bond.add(other_bonds[0].id_)
            self._seen_atom.add(next_atom.id_)
            self.number_of_double_bonds_in_series += 1
            return self.check_for_multiple_double_bonds(next_atom, other_bonds[0])

        return atom

    def determine_up_down_atom_left(self, start_atom: Atom):
        bond_ez = get_bond_with_ez_symbol(start_atom)
        if bond_ez is None:
            self.setup_successful = False
            return  # no stereo bonds '/' or '\\'

        atom_ez = get_other_node(bond_ez, start_atom)
        self._seen_bond.add(bond_ez.id_)
        self._seen_atom.add(atom_ez.id_)

        # set up path
        path = TraversalPath(atom_ez, bond_ez)
        if bond_ez.symbol == "/":
            if atom_ez.parent is start_atom.parent:
                self.up = path  # main chain
            else:
                self.down = path  # branch
        else:  # "\\"
            if atom_ez.parent is start_atom.parent:
                self.down = path  # main chain
            else:
                self.up = path  # branch

        self.figure_out_other_bond(start_atom, bond_ez)

    def determine_up_down_atom_right(self, start_atom: Atom):
        bond_ez = get_bond_with_ez_symbol(start_atom)
        if bond_ez is None:
            self.setup_successful = False
            return  # no stereo bonds '/' or '\\'

        atom_ez = get_other_node(bond_ez, start_atom)
        self._seen_bond.add(bond_ez.id_)
        self._seen_atom.add(atom_ez.id_)

        # set first atom that has stereo symbol
        path = TraversalPath(atom_ez, bond_ez)
        if bond_ez.symbol == "/":
            self.down = path
        else:  # "\\"
            self.up = path

        self.figure_out_other_bond(start_atom, bond_ez)

    def figure_out_other_bond(self, start_atom: Atom, bond_ez: Bond):
        # figure out other bond
        other_bond = get_other_bonds(start_atom, (self.double_bond, bond_ez))
        if len(other_bond) == 1:
            other_atom = get_other_node(other_bond[0], start_atom)
            if self.down is None:
                self.down = TraversalPath(other_atom, other_bond[0])
            else:
                self.up = TraversalPath(other_atom, other_bond[0])

            self.check_for_winner()

        elif len(other_bond) == 0:
            # other bond must be implicit hydrogens or terminal
            if self.down is None:
                self.winner = TraversalOutcomes.up
            else:
                self.winner = TraversalOutcomes.down

        else:
            raise NotImplementedError("more than 2 bonds coming off double bond atom")

    def step(self):
        # rule 1 higher atomic mass has higher priority
        # rule 2 propagate down chain till there is a difference

        # traverse the graph one step
        self.up.next_step()
        self.down.next_step()

        # see if there is a winner
        self.check_for_winner()

    def check_for_winner(self):
        score_up = get_score(self.up.current_node)
        score_down = get_score(self.down.current_node)
        if score_up == score_down:
            if score_up is None:
                self.winner = TraversalOutcomes.tie
            return

        if score_up > score_down:
            self.winner = TraversalOutcomes.up
        else:
            self.winner = TraversalOutcomes.down


def check_setup(left, right) -> bool:
    """
    Check to see if there were reasons to terminate early from traversal setup

    Parameters
    ----------
    left:

    right:


    Returns
    -------
    result:
        True: terminate early, False: continue with graph traversal

    """
    # check for 'odd' number of double bonds in a row
    if (left.number_of_double_bonds_in_series + right.number_of_double_bonds_in_series + 1) % 2 == 0:
        # left + right + 1 (starting double bond) = even  --> return None for E/Z
        logging.warning("Double bond stereo-chemistry notation '/' and '\\' should "
                        f"not be used on even number of double bonds; they are axial chiral and use @/@@ notation on "
                        f"center atom.\nBigSMILES:{left.double_bond.parent}")
        return True

    return False


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

    * Depth-first-search (DFS) prioritizing heavy atoms and following multiple branches if atom count equal.
    Starts out at the double bond and heads left and right,
    then will branch into up and down for a total of 4 simultaneous DFS. At each step the higher atomic number
    element is chosen (if same, backtracking is noted)
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

    if not left_side.setup_successful or not right_side.setup_successful or check_setup(left_side, right_side):
        return None

    logging.debug(f" Starting E/Z graph traversal: {double_bond.root.to_string(show_repr=(Atom,))}"
                  f"\n-------------------------------")
    logging.debug(" step, left_up, left_down, right_up, right_down")
    logging.debug(f" {0:4}, {repr(left_side.up.current_node):>7}, {repr(left_side.down.current_node):>9}, "
                  f"{repr(right_side.up.current_node):>8}, {repr(right_side.down.current_node):>10}")

    # perform graph traversal checking each step to see if done
    for i in range(200):  # 200 is the max depth of traversal
        if not left_side.done:
            left_side.step()

        if not right_side.done:
            right_side.step()

        logging.debug(f" {i:4}, {repr(left_side.up.current_node):>7}, {repr(left_side.down.current_node):>9}, "
                      f"{repr(right_side.up.current_node):>8}, {repr(right_side.down.current_node):>10}")
        if left_side.done and right_side.done:
            logging.debug("Traversal Done \n")
            break  # traversal done

    else:
        logging.warning(f'Max traversal depth reached on E/Z determination of: '
                        f'\nBigSMILES: {double_bond.parent} '
                        f'\nBond: {repr(double_bond)}')
        return None

    # process results
    if left_side.winner is TraversalOutcomes.tie and right_side.winner is TraversalOutcomes.tie:
        logging.warning("Double bond stereo-chemistry notation '/' and '\\' should "
                        "not be used on symmetric molecules.")
        return None
    if (left_side.winner is TraversalOutcomes.tie or right_side.winner is TraversalOutcomes.tie) or \
            (left_side.winner == right_side.winner):
        return "Z"

    return "E"
