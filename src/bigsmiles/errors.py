

class BigSMILESError(Exception):
    def __init__(self, text: str, error=None):
        if isinstance(error, BigSMILESError):
            text += '\n\t'
            text += error.text.replace("\t", "\t\t")
        self.text = text

    def __str__(self):
        return self.text


class TokenizeError(BigSMILESError):
    """ Raised when an error occurs tokenizing a BigSMILES string. """


class ConstructorError(BigSMILESError):
    """ Raised when an error occurs during the construction of a BigSMILES string. """


class ValidationError(BigSMILESError):
    """ Raised when validating a SMILES/BigSMILES string is syntactically correct. """


ERRORS = (BigSMILESError, TokenizeError)
