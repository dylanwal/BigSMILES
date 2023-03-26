from __future__ import annotations


class BigSMILESError(Exception):
    """ base BigSMILES error """

    def __init__(self, text: str, error=None):
        if isinstance(error, BigSMILESError):
            text += '\n\t'
            text += error.text.replace("\t", "\t\t")
        self.text = text

    def __str__(self):
        return self.text


class TokenizeError(BigSMILESError):
    """ raised when an error occurs tokenizing a BigSMILES string. """


class ConstructorError(BigSMILESError):
    """ raised when an error occurs during the construction of a BigSMILES string. """


class ValidationError(BigSMILESError):
    """ raised when validating a SMILES/BigSMILES string is syntactically correct. """


class DistributionError(BigSMILESError):
    """ raise when issue with MW distribution """


class MolecularFormulaError(BigSMILESError):
    """ raise when issue with calculating molecular formula """
