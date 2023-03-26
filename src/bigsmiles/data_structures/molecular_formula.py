from __future__ import annotations

import re

import bigsmiles.errors as errors
import bigsmiles.reference_data.chemical_data as chemical_data

# TODO: take allow isotopes in calculation - get from BigSMILES


class MolecularFormula:
    """
    A list of the number of atoms present in the molecule or molecular fragment.

    """
    __slots__ = ("_elements", "_formula", "_elemental_analysis", "_molar_mass")

    def __init__(self, element_dict: dict[str, int] = None):
        self._elements = element_dict
        self._formula: str | None = None
        self._elemental_analysis: dict[str, float] | None = None
        self._molar_mass: int | float | None = None

    def __str__(self):
        return self.formula

    def __repr__(self):
        return self.formula

    def __getitem__(self, key: str):
        return self._elements[key]

    def __iter__(self) -> tuple[str, int]:
        for element, num in self._elements.items():
            yield element, num

    @property
    def formula(self) -> str:
        """
        string representation of molecular formula (Hills system)

        !!! warning

            '{}' is appended to the end when a stochastic object is present

        Notes
        -----
        Hills system (Hill notation) is the default ordering where cabon followed by hydrogen is first then
        all other elements in alphabetical order. If no carbon is present, all elements (including hydrogen) are listed
        alphabetically. Also, single-letter elements coming before two-letter symbols.
        Reference: J. Am. Chem. Soc. 1900, 22, 8, 478. https://pubs.acs.org/doi/abs/10.1021/ja02046a005
        """
        if self._formula is None:
            self._formula = elements_to_formula_string(self._elements)

        return self._formula

    @property
    def elements(self) -> dict[str, int]:
        """ dictionary representation of molecular formula;  {symbol: count} """
        return self._elements

    @property
    def contains_stochastic_object(self) -> bool:
        """ True: stochastic object within molecular formula """
        if "{}" in self._elements:
            return True
        return False

    @property
    def molar_mass(self) -> int | float:
        """
        mass for one mole of molecule or molecular fragment

        Returns
        -------
        molar_mass:
            molar mass in g/mol

        !!! warning
            Will return a number with stochastic object inside. Use 'contains_stochastic_object' to check.

        """
        if self._molar_mass is None:
            self._elemental_analysis, self._molar_mass = compute_element_analysis_and_molar_mass(self.elements)

        return self._molar_mass

    @property
    def elemental_analysis(self) -> dict[str, int | float]:
        """
        mass fraction of each element; {symbol: mass fraction}

        Returns
        -------
        elemental_analysis:
            mass fraction of each element

        !!! warning
            Will return a number with stochastic object inside. Use 'contains_stochastic_object' to check.
        """
        if self._elemental_analysis is None:
            self._elemental_analysis, self._molar_mass = compute_element_analysis_and_molar_mass(self.elements)

        return self._elemental_analysis


def compute_element_analysis_and_molar_mass(elements: dict[str, int]) -> tuple[dict[str, float], int | float]:
    element_analysis = dict()
    for element in elements:
        try:  # try loop used for speed benefit over if statement
            element_analysis[element] = chemical_data.symbol_to_atomic_mass[element] * elements[element]
        except KeyError as e:
            if element == "{}":
                continue
            raise e

    molar_mass = sum([v for v in element_analysis.values()])

    for k in element_analysis:
        element_analysis[k] = element_analysis[k] / molar_mass

    return element_analysis, molar_mass


def parse_molecular_formula_strings(formula: str) -> MolecularFormula:
    """
    Parse molecular formula string to MolecularFormula

    Parameters
    ----------
    formula:
        molecular formula string (only small molecules)

    Returns
    -------
    molecular_formula:
        molecular_formula object

    Examples
    --------
    >>> parse_molecular_formula_strings("CCOCC")

    """
    split_formula = split_formula_string(formula)
    elements = split_formula_to_element_dict(split_formula)
    return MolecularFormula(elements)


def split_formula_string(formula: str) -> list[str]:
    """
    Split formula up into letters and numbers.

    Parameters
    ----------
    formula:
        string formula

    Returns
    -------
    split_formula:

    Examples
    --------
    >>> split_formula_string("C1H3OH63Cr2CCCOOO")
    ['C', '1', 'H', '3', 'O', 'H', '63', 'Cr', '2', 'C', 'C', 'C', 'O', 'O', 'O']

    """
    pattern = r'[A-Z]{1}[a-z]*|\d+'
    return re.findall(pattern, formula)


def split_formula_to_element_dict(split_formula: list[str]) -> dict[str, int]:
    """
    Given a split chemical formula return a dictionary of element_symbols.

    Parameters
    ----------
    split_formula:
        split chemical formula

    Returns
    -------
    elements:
        {"symbol": "count"}

    Examples
    --------
    >>> split_formula_to_element_dict(['O', 'C', 'O'])
    {"C": 1, "O": 2}
    >>> split_formula_to_element_dict(['C', '2', 'O', 'N', '4', 'F'])
    {"C": 2, "O": 1, "N": 4, "F": 1}

    """
    # put data into dictionary
    elements = {}
    for index, entry in enumerate(split_formula):
        if not entry.isnumeric():
            if entry not in chemical_data.element_symbols:
                raise errors.MolecularFormulaError(f"Invalid symbol in chemical formula. Invalid: {entry}")

            if entry in elements:
                if index + 1 < len(split_formula):
                    num_ = split_formula[index + 1]
                else:
                    num_ = "1"
                if num_.isnumeric():
                    elements[entry] = int(num_) + elements[entry]
                else:
                    elements[entry] = 1 + elements[entry]
            else:
                if index + 1 < len(split_formula):
                    num_ = split_formula[index + 1]
                else:
                    num_ = "1"
                if num_.isnumeric():
                    elements[entry] = int(num_)
                else:
                    elements[entry] = 1

    return elements


def elements_to_formula_string(elements: dict[str, int], order: list[str] = None) -> str:
    """
    Converts symbol dictionary into a reduced string.

    Parameters
    ----------
    elements
    order

    Returns
    -------
    chemical_formula:
        chemical formula

    Examples
    --------
    >>> elements_to_formula_string({"C": 2, "O": 1, "N": 4, "F": 1})
    "C2ONF"

    Notes
    -----
    Hills system (Hill notation) is the default ordering where cabon followed by hydrogen is first then
    all other elements in alphabetical order. If no carbon is present, all elements (including hydrogen) are listed
    alphabetically. Also, single-letter elements coming before two-letter symbols.
    Reference: J. Am. Chem. Soc. 1900, 22, 8, 478. https://pubs.acs.org/doi/abs/10.1021/ja02046a005

    """
    # Set order
    if order is None:
        # default is hill system
        if "C" in elements.keys():
            order = set_element_order_in_formula(["C", "H"])
        else:
            order = set_element_order_in_formula()
    else:
        order = set_element_order_in_formula(order)

    # Generate chemical formula
    chemical_formula = ""
    for element in order:
        if element in elements.keys() and (elements[element] is not None or elements[element] == 0):
            if elements[element] == 1:
                chemical_formula += element
            else:
                chemical_formula += element + str(elements[element])

    if "{}" in elements:
        chemical_formula += "{}" + str(elements["{}"] if elements["{}"] is not 1 else "")

    return chemical_formula


def set_element_order_in_formula(order: list[str] = None) -> list[str]:
    """Given a chemical order of element_symbols, validate it/ and complete it. Default is alphabetically."""
    base = list(chemical_data.elements_ordered)
    if order is None:
        return base

    for k in reversed(order):
        if k in base:
            base.remove(k)
            base.insert(0, k)
        else:
            raise errors.MolecularFormulaError(f"Invalid symbol in 'order' for chemical formula. Invalid: {k}")

    return base
