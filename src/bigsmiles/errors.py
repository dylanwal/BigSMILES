

class BigSMILESError(Exception):
    def __init__(self, text: str, error=None):
        if isinstance(error, BigSMILESError):
            text += '\n\t'
            text += error.text.replace("\t", "\t\t")
        self.text = text

    def __str__(self):
        return self.text

