import math

import pytest

import bigsmiles.errors as errors
import bigsmiles.reference_data.chemical_data as chemical_data
import bigsmiles.data_structures.molecular_formula as molecular_formula
import bigsmiles


def test_set_element_order_in_formula():
    order = molecular_formula.set_element_order_in_formula()
    for i, ii in zip(order, chemical_data.elements_ordered):
        assert i == ii


def test_set_element_order_in_formula2():
    specific_order = ["C", "H"]
    order = molecular_formula.set_element_order_in_formula(specific_order)
    answer = list(chemical_data.elements_ordered)
    answer.remove("C")
    answer.remove("H")
    answer = specific_order + answer
    for i, ii in zip(order, answer):
        assert i == ii


def test_set_element_order_in_formula_invalid_order():
    order = ["H", "O", "Zn", "PP"]  # 'PP' invalid
    with pytest.raises(errors.MolecularFormulaError) as _:
        molecular_formula.set_element_order_in_formula(order)


def test_elements_to_formula_string():
    result = molecular_formula.elements_to_formula_string({"C": 2, "O": 1, "N": 4, "F": 1, "H": 4})
    assert result == "C2H4FN4O"


def test_split_formula_to_element_dict():
    result = molecular_formula.split_formula_to_element_dict(['C', '2', 'O', 'N', '4', 'F'])
    answer = {"C": 2, "O": 1, "N": 4, "F": 1}
    for k in result:
        assert result[k] == answer[k]


def test_split_formula_string():
    result = molecular_formula.split_formula_string("C1H3OH63Cr2CCCOOO")
    answer = ['C', '1', 'H', '3', 'O', 'H', '63', 'Cr', '2', 'C', 'C', 'C', 'O', 'O', 'O']
    for i, ii in zip(result, answer):
        assert i == ii


cases_smiles = [
    ["CCCCCC",
     {
         "mw": 86.18,
         "formula": "C6H14",
         "elemental_analysis": {"C": 0.8362, "H": 0.1638}
     }
     ],
    ["c1ccccc1",
     {
         "mw": 78.11,
         "formula": "C6H6",
         "elemental_analysis": {"C": 0.9226, "H": 0.0774}
     }
     ],
    [r"F/C(C)=C(C(OC(NBr)C)C(C)(F)C)\C(OC(N)C)C(F)=C",
     {
         "mw": 389.26,
         "formula": "C14H24BrF3N2O2",
         "elemental_analysis": {"C": 0.4320, "H": 0.0621, "Br": 0.2053, "F": 0.1464, "N": 0.0720, "O": 0.0822}
     }
     ],
    ["c1ccccc1",
     {
         "mw": 78.11,
         "formula": "C6H6",
         "elemental_analysis": {"C": 0.9226, "H": 0.0774}
     }
     ]
]


@pytest.mark.parametrize("case", cases_smiles)
def test_smiles(case: list[str, dict]):
    bigsmiles_string, answers = case
    mol = bigsmiles.BigSMILES(bigsmiles_string)

    assert mol.molecular_formula.formula == answers["formula"]
    assert math.isclose(mol.molar_mass, answers["mw"], rel_tol=1e-3)
    elemental_analysis = answers["elemental_analysis"]
    for k in elemental_analysis:
        assert math.isclose(mol.molecular_formula.elemental_analysis[k], elemental_analysis[k], rel_tol=1e-3)


def test_bigsmiles():
    mol = bigsmiles.BigSMILES("[H]{[$][$]CC[$][$]}[H]")

    assert mol.molecular_formula.contains_stochastic_object is True
    assert mol.molecular_formula.formula == "H2{}"
    assert math.isclose(mol.molar_mass, 2.016, rel_tol=1e-3)
    elemental_analysis = {"H": 1}
    for k in elemental_analysis:
        assert math.isclose(mol.molecular_formula.elemental_analysis[k], elemental_analysis[k], rel_tol=1e-3)

    stoch_obj = mol.nodes[2]
    stoch_frag = stoch_obj.nodes[0]

    assert stoch_frag.molecular_formula.formula == "C2H4"
    assert math.isclose(stoch_frag.molar_mass, 28.0540, rel_tol=1e-3)


