"""

Validation on SMILES or BigSMILES string.

"""
import re

from bigsmiles.errors import ValidationError


def run_string_validation(text: str) -> str:
    """ Main entry point for pre-construction validation. """
    branch_symbol_validation(text)
    stochastic_object_validation(text)
    brackets_validation(text)
    # text = validate_ring_numbering(text)  # Faster to do after tokenization
    return text


def branch_symbol_validation(text: str):
    """ Confirm that there is an equal number of '(' to ')'. """
    count = text.count('(') - text.count(')')
    if count != 0:
        if count > 0:
            raise ValidationError(f"Invalid BigSMILES. Missing {count} closing branch symbols ')'. "
                                  f"\n Invalid string: {text}")
        raise ValidationError(f"Invalid BigSMILES. Missing {count} opening branch symbols '('. "
                              f"\n Invalid string: {text}")


def stochastic_object_validation(text: str):
    """ Confirm that there is an equal number of '{' to '}'. """
    count = text.count('{') - text.count('}')
    if count != 0:
        if count > 0:
            raise ValidationError(f"Invalid BigSMILES. Missing {count} closing stochastic object symbols" + " '}'. "
                                                                                                            f"\n Invalid string: {text}")
        raise ValidationError(f"Invalid BigSMILES. Missing {count} opening stochastic object symbols " + "'{'. "
                                                                                                         f"\n Invalid string: {text}")


def brackets_validation(text: str):
    """ Confirm that there is an equal number of '[' to ']'. """
    count = text.count('[') - text.count(']')
    if count != 0:
        if count > 0:
            raise ValidationError(f"Invalid BigSMILES. Missing {count} closing bracket symbols ']'. "
                                  f"\n Invalid string: {text}")
        raise ValidationError(f"Invalid BigSMILES. Missing {count} opening bracket symbols '['. "
                              f"\n Invalid string: {text}")


ring_digit_pattern = re.compile(r"(?<![H+><$-]|\[)(%\d{2}|\d)")
extract_bracket_pattern = re.compile(r"(\[)(.*?)(\])")

ascii_lowercase = 'abcdefghijklmnopqrstuvwxyz'
len_ascii_lowercase = len(ascii_lowercase)


def get_next_letter(num: int) -> str:
    num_sym = int(num/len_ascii_lowercase)
    mod_sym = num % len_ascii_lowercase

    if num == 0:
        sym = ascii_lowercase[mod_sym]
    else:
        sym = ascii_lowercase[mod_sym] * (num_sym + 1)

    return "<" + sym + ">"


def validate_ring_numbering(text: str, renumber: bool = True):
    """
    Confirm that there is only two of a number.

    Note: Faster to do after tokenization
    """
    # No rings
    if text.count("1") == 0:
        return text

    # remove atoms with special stuff as it mess with ring parsing
    extract = dict()
    split_list = re.split(extract_bracket_pattern, text)
    for i in range(len(split_list)):
        if split_list[i] == "[":
            replacement_symbol = get_next_letter(i)
            extract[replacement_symbol] = split_list[i + 1]
            split_list[i + 1] = replacement_symbol

    text = "".join(split_list)

    ring_numbers = [int(num.replace("%", "")) for num in ring_digit_pattern.findall(text)]
    ring_number_set = set(ring_numbers)

    # check for more than two
    for num in list(ring_number_set):
        count = ring_numbers.count(num)
        if count != 2:
            if count == 1:
                raise ValidationError(f"Invalid BigSMILES. Missing ring index. Only one found for: '{num}'. "
                                      f"\n Invalid string: {text}")
            elif renumber and count % 2 == 0:
                # renumber duplicates
                text = _renumber_rings(text, num, max(ring_number_set) + 1)
                for i in range(1, int(count / 2) + 1):
                    ring_number_set.add(max(ring_number_set) + i)

            else:
                raise ValidationError(f"Invalid BigSMILES. More than 2 ring index found for: '{num}'. "
                                      f"\n Invalid string: {text}")

    # do replacement
    for k, v in extract.items():
        text = text.replace(k, v)
    return text


def _renumber_rings(text: str, num: int, new_num: int) -> str:
    if num > 9:
        pattern = r"(?<![H+><$-]|\[)" + f"(%{num})"
    else:
        pattern = r"(?<![H+><$-%]|\[)" + f"({num})"

    hit_index = [(m.start(0), m.end(0)) for m in re.finditer(pattern, text)]

    new_num += int(len(hit_index) / 2 - 2)
    first = True
    for hit in reversed(hit_index[2:]):  # skip the first two
        if new_num > 9:
            text = text[:hit[0]] + "%" + str(new_num) + text[hit[1]:]
        else:
            text = text[:hit[0]] + str(new_num) + text[hit[1]:]

        if first:
            first = False
        else:
            first = True
            new_num -= 1

    return text
