"""

Rules for determining double bond E/Z stereo-chemistry

    * rule 1: higher atomic mass has higher priority
    * rule 2: propagate down chain till there is a difference

Approach starting at the double bond perform a graph traversal in all 4 direction simultaneously.

Note: stochastic objects are scored as 150 and treated as a single atom.

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
    for bond in atom.bonds:
        if bond.symbol == "/" or bond.symbol == "\\":
            return bond

    return None


def get_score(node: Atom | StochasticObject | BondDescriptorAtom | None) -> int:
    if node is None:
        return 0
    elif isinstance(node, Atom):
        return chemical_data.symbol_to_atomic_number[node.element]

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
            if self.current_node.bond_left == self.prior_edge:
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

    def __init__(self, double_bond: Bond, right_left: str):
        self.double_bond = double_bond
        self.right_left = right_left

        self.up: TraversalPath | None = None
        self.down: TraversalPath | None = None

        self.winner: TraversalOutcomes | None = None

        if right_left == "left":
            self.determine_up_down_atom_left()
        else:
            self.determine_up_down_atom_right()

    @property
    def done(self) -> bool:
        if self.winner is None:
            return False
        return True

    def determine_up_down_atom_left(self):
        atom_from_double_bond = self.double_bond.atom1
        bond_ez = get_bond_with_ez_symbol(atom_from_double_bond)
        atom_ez = get_other_node(bond_ez, atom_from_double_bond)

        # set first atom that has stereo symbol
        path = TraversalPath(atom_ez, bond_ez)
        if bond_ez.symbol == "/":
            if atom_ez.parent is atom_from_double_bond.parent:
                self.up = path  # main chain
            else:
                self.down = path  # branch
        else:  # "\\"
            if atom_ez.parent is atom_from_double_bond.parent:
                self.down = path  # main chain
            else:
                self.up = path  # branch

        self.figure_out_other_bond(atom_from_double_bond, bond_ez)

    def determine_up_down_atom_right(self):
        atom_from_double_bond = self.double_bond.atom2
        bond_ez = get_bond_with_ez_symbol(atom_from_double_bond)
        atom_ez = get_other_node(bond_ez, atom_from_double_bond)

        # set first atom that has stereo symbol
        path = TraversalPath(atom_ez, bond_ez)
        if bond_ez.symbol == "/":
            if atom_ez.parent is atom_from_double_bond.parent:
                self.down = path  # main chain
            else:
                self.up = path  # branch
        else:  # "\\"
            if atom_ez.parent is atom_from_double_bond.parent:
                self.up = path  # main chain
            else:
                self.down = path  # branch

        self.figure_out_other_bond(atom_from_double_bond, bond_ez)

    def figure_out_other_bond(self, atom_from_double_bond: Atom, bond_ez: Bond):
        # figure out other bond
        other_bond = get_other_bonds(atom_from_double_bond, (self.double_bond, bond_ez))
        if other_bond:
            other_atom = get_other_node(other_bond[0], atom_from_double_bond)
            if self.down is None:
                self.down = other_atom
            else:
                self.up = other_atom

            self.check_for_winner()
            return

        # other bond must be to [H]
        if self.down is None:
            self.winner = TraversalOutcomes.up
        else:
            self.winner = TraversalOutcomes.down

    def step(self):
        # rule 1 higher atomic mass has higher priority
        # rule 2 propagate down chain till there is a difference

        # traverse the graph one step
        self.up.next_step()
        self.down.next_step()

        # see if there is a winner
        self.check_for_winner()

    def check_for_winner(self):
        score_up = get_score(self.up)
        score_down = get_score(self.down)
        if score_up == score_down:
            if score_up is None:
                self.winner = TraversalOutcomes.tie
            return
        if score_up > score_down:
            self.winner = TraversalOutcomes.up
        self.winner = TraversalOutcomes.down


def get_double_bond_ez(double_bond: Bond) -> str | None:
    """
    graph traversal in 4 bond coming off double bond

    """
    left_side = TraversalSide(double_bond, "left")
    right_side = TraversalSide(double_bond, "right")

    for i in range(100):  # 100 is the max depth of traversal
        if not left_side.done:
            left_side.step()

        if not right_side:
            right_side.step()

        if left_side.done and right_side.done:
            break

    if left_side.winner is TraversalOutcomes.tie and right_side.winner is TraversalOutcomes.tie:
        logging.warning("Double bond stereo-chemistry notation '/' and '\\' should "
                        "not be used on symmetric molecules.")
        return None

    if left_side.winner is TraversalOutcomes.tie or right_side.winner is TraversalOutcomes.tie:
        return "Z"

    if left_side.winner == right_side.winner:
        return "Z"

    return "E"
