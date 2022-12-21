

class BigSMILESError(Exception):
    def __init__(self, text: str, optional: str = None):
        if optional is not None:
            text += "(" + optional + ")"
        super().__init__(text)

