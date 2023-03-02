
import logging

from bigsmiles.errors import ValidationError
from bigsmiles.constructors.tokenizer import Token, TokenKind


def run_token_validation(tokens: list[Token]) -> list[Token]:
    tokens = validate_ring_numbering(tokens)

    return tokens


def validate_ring_numbering(tokens: list[Token], renumber: bool = True) -> list[Token]:
    """ Confirm that there is only two of a number. """

    # get ring_numbers
    ring_numbers = dict()  # key: ring number; value: list[index in tokens]
    for i, token in enumerate(tokens):
        if token.kind is TokenKind.Ring:
            value = int(token.value)
        elif token.kind is TokenKind.Ring2:
            value = int(token.value.replace("%", ""))
        else:
            continue

        if value in ring_numbers:
            ring_numbers[value].append(i)
        else:
            ring_numbers[value] = [i]

    # No rings
    if not ring_numbers:
        return tokens

    # check for more than two
    max_index = max([k for k in ring_numbers.keys()]) + 1
    for key in list(ring_numbers.keys()):
        count = len(ring_numbers[key])
        if count == 2:
            continue
        elif count == 1:
            raise ValidationError(f"Invalid BigSMILES. Missing ring index. Only one found for: '{key}'. ")
        elif renumber and count % 2 == 0:
            logging.warning(f"Duplicate ring index detected ({key}) and fixed for: "
                            f"{''.join([str(t.value) for t in tokens])}")
            # renumber duplicates
            for i, index_ in enumerate(ring_numbers[key]):
                if int(i/2) == 0:
                    continue  # don't change the first set

                tokens[index_].value = str(max_index) if max_index < 10 else "%" + str(max_index)
                if i % 2 != 0:
                    max_index += 1
        else:
            raise ValidationError(f"Invalid BigSMILES. More than 2 ring index found for: '{key}'. ")

    return tokens

