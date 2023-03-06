from __future__ import annotations

import bigsmiles.errors as errors
import bigsmiles.reference_data.chemical_data as chemical_data
from bigsmiles.data_structures.bigsmiles import Bond, Atom




def get_other_atom(bond: Bond, not_this_atom: Atom) -> Atom | None:
    for atom in bond:
        if atom is not not_this_atom:
            return atom

    return None


def get_other_bonds(atom: Atom, not_this_bond: tuple[Bond, ...] | list[Bond, ...]) -> list[Bond, ...]:
    bonds = []
    for bond in atom.bonds:
        if bond not in not_this_bond:
            bonds.append(bonds)

    return bonds


def get_bond_with_ez_symbol(atom: Atom) -> Bond | None:
    for bond in atom.bonds:
        if bond.symbol == "/" or bond.symbol == "\\":
            return bond

    return None


class DoubleBondSide:

    def __init__(self, double_bond: Bond, right_left: str):
        self.double_bond = double_bond
        self.right_left = right_left

        self.up: Atom | None = None
        self.down: Atom | None = None

        self.winner = None
        self.done = False

        if right_left == "left":
            self.determine_up_down_atom_left()
        else:
            self.determine_up_down_atom_right()

    def __bool__(self):
        return self.done

    def determine_up_down_atom_left(self):
        atom_from_double_bond = self.double_bond.atom1
        bond_ez = get_bond_with_ez_symbol(atom_from_double_bond)
        atom_ez = get_other_atom(bond_ez, atom_from_double_bond)

        # set first atom that has stereo symbol
        if bond_ez.symbol == "/":
            if atom_ez.parent is atom_from_double_bond.parent:
                self.up = atom_ez  # main chain
            else:
                self.down = atom_ez  # branch
        # "\\"
        if atom_ez.parent is atom_from_double_bond.parent:
            self.down = atom_ez  # main chain
        else:
            self.up = atom_ez  # branch

        self.figure_out_other_bond(atom_from_double_bond, bond_ez)

    def determine_up_down_atom_right(self):
        atom_from_double_bond = self.double_bond.atom2
        bond_ez = get_bond_with_ez_symbol(atom_from_double_bond)
        atom_ez = get_other_atom(bond_ez, atom_from_double_bond)

        # set first atom that has stereo symbol
        if bond_ez.symbol == "/":
            if atom_ez.parent is atom_from_double_bond.parent:
                self.down = atom_ez  # main chain
            else:
                self.up = atom_ez  # branch
        # "\\"
        if atom_ez.parent is atom_from_double_bond.parent:
            self.up = atom_ez  # main chain
        else:
            self.down = atom_ez  # branch

        self.figure_out_other_bond(atom_from_double_bond, bond_ez)

    def figure_out_other_bond(self, atom_from_double_bond: Atom, bond_ez: Bond):
        # figure out other bond
        other_bond = get_other_bonds(atom_from_double_bond, (self.double_bond, bond_ez))
        if other_bond:
            other_atom = get_other_atom(other_bond[0], atom_from_double_bond)
            if self.down is not None:
                self.down = other_atom
            else:
                self.up = other_atom

            return

        # other bond must be to [H]
        if self.down is not None:
            self.winner = self.up
        else:
            self.winner = self.down

        self.done = True

    def step(self):
        # rule 1 higher atomic mass has higher priority
        # rule 2 probagate down chain till there is a difference
        if self.right_left == "left":
            self.double_bond.atom1


def get_double_bond_ez(double_bond: Bond) -> str:
    """
    graph traversal in 4 bond coming off double bond

    """
    left_side = DoubleBondSide(double_bond, "left")
    right_side = DoubleBondSide(double_bond, "right")

    for i in range(100):  # 100 is the max depth of traversal
        if left_side:
            left_side.step()

        if right_side:
            right_side.step()

        if not left_side and not right_side:
            break

    if left_side.winner == right_side.winner:
        return "z"

    return "e"
