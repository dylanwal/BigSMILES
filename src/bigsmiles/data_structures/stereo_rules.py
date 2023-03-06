from __future__ import annotations

import bigsmiles.errors as errors
import bigsmiles.reference_data.chemical_data as chemical_data
from bigsmiles.data_structures.bigsmiles import Bond, Atom

# rule 1 higher atomic mass has higher priority
# rule 2 probagate down chain till there is a difference


class DoubleBondSide:

    def __init__(self, double_bond: Bond, right_left: str):
        self.double_bond = double_bond
        self.right_left = right_left

        self.up = None
        self.down = None

        self.winner = None
        self.done = False

    def __bool__(self):
        return self.done

    def determine_up_down(self):
        if self.right_left == "left":
            first_bond = self.double_bond.atom1.bonds[0]  # main chain
            second_bond = self.double_bond.atom1.bonds[1] if self.double_bond.atom1.bonds[1] is not self.double_bond else None
        else:
            first_bond = double_bond.atom1.bonds[0]  # branch of main chain
            second_bond = double_bond.atom1.bonds[1] if double_bond.atom1.bonds[1] is not double_bond else None


            if second is None:
                self.main = first
                self.branch = None
            else:
                self.main = second
                self.branch = first
        else:
            raise ValueError("Invalid value for 'right_left'. ")

    def step(self):
        if self.right_left == "left":
            self.double_bond.atom1


def get_double_bond_ez(double_bond: Bond):
    """
    graph traversal in 4 bond coming off double bond

    """
    return None
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
