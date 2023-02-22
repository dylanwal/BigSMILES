import warnings


def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'

warnings.formatwarning = custom_formatwarning

from bigsmiles.validation import pre_validation
from bigsmiles.tokenizer import Token, TokenKind, tokenize
from bigsmiles.constructor import States
from bigsmiles.constructor_str import BigSMILESStringConstructor
from bigsmiles.errors import BigSMILESError


def map_atom(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token):
    if constructor.state is States.start:
        constructor.add_atom(token.value)

    elif constructor.state in (States.atom, States.bond_descriptor, States.stochastic_object_end,
                               States.branch_start, States.ring):
        constructor.add_bond_atom_pair("", token.value)  # add single bond
    elif constructor.state in (States.stochastic_fragment,):
        constructor.add_atom(token.value)

    else:
        raise BigSMILESError("Something went wrong in 'map_atom'.")  # should not be hit


def map_bond(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token):
    try:
        next_token = tokens.pop(0)
    except IndexError:
        raise BigSMILESError("Bond can't be at the end of a BigSMILES string.")

    if next_token.kind in (TokenKind.Atom, TokenKind.AtomExtend):
        constructor.add_bond_atom_pair(token.value, next_token.value)
    elif next_token.kind is not TokenKind.BondDescriptor:
        constructor.add_bond_bonding_descriptor_pair_str(token.value, next_token.value)
    elif next_token.kind is not TokenKind.StochasticStart:
        map_stochastic_object_start(constructor, tokens, next_token, token)
    else:
        raise BigSMILESError(f"Bonds must be followed by an Atom, Bonding Descriptor, or Stocastic Object. "
                             f"Parse completed: {constructor.bigsmiles}, issue token: {token}")


def map_bond_descriptor(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token):
    if tokens[0].kind is TokenKind.StochasticEnd:
        constructor.close_stochastic_fragment()
        constructor.close_stochastic_object_str(token.value)
        tokens.pop(0)
        return

    if constructor.state is States.stochastic_fragment:
        # first StochasticFragment symbol
        constructor.add_bonding_descriptor_str(token.value)
        return

    if constructor.state is States.branch_start:
        # first Branch symbol
        # check to make sure the branch closes immediately
        if tokens[0].kind != TokenKind.BranchEnd:
            raise BigSMILESError("If Bonding Descriptors is the first Branch symbol, it can only be followed by "
                                 "Branch End.")
    constructor.add_bond_bonding_descriptor_pair_str("", token.value)


def map_branch_start(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token):
    constructor.open_branch()


def map_branch_end(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token):
    constructor.close_branch()


def map_ring(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token):
    if constructor.state is not States.atom:
        raise BigSMILESError(f"Ring number must follow atoms.")

    constructor.add_ring(int(token.value.replace('%', '')))


def map_stochastic_object_start(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token,
                                bond_token: Token = None):
    try:
        next_token = tokens.pop(0)
    except IndexError:
        raise BigSMILESError("Stochastic objects must begin with a bond descriptor (or implicit bonding descriptor).")

    if next_token.kind not in (TokenKind.ImplictEndGroup, TokenKind.BondDescriptor):
        raise BigSMILESError(f"Stochastic object starts must be followed an explict or implicit end group.")

    if constructor.state is States.start:
        constructor.open_stochastic_object_str(next_token.value)
    else:
        if bond_token is None:
            constructor.open_stochastic_object_with_bond_str("", next_token.value)
        else:
            constructor.open_stochastic_object_with_bond_str(bond_token.value, next_token.value)


def map_stochastic_object_end(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token) -> TokenKind:
    raise BigSMILESError("Stochastic objects should end with bonding descriptor (or implicit bonding description)")


def map_bond_seperator(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token):
    constructor.close_open_stochastic_fragment()


def NotImplementedFunc(*args, **kwargs):
    raise NotImplementedError


def SkipSymbol(constructor: BigSMILESStringConstructor, tokens: list[Token], token: Token):
    warnings.warn(f"Symbol skipped: {token.value}")


map_tokens = {
    TokenKind.Bond: map_bond,
    TokenKind.Atom: map_atom,
    TokenKind.Aromatic: map_atom,
    TokenKind.AtomExtend: map_atom,
    TokenKind.BranchStart: map_branch_start,
    TokenKind.BranchEnd: map_branch_end,
    TokenKind.Ring: map_ring,
    TokenKind.Ring2: map_ring,
    TokenKind.BondEZ: SkipSymbol,
    TokenKind.Mix: NotImplementedFunc,
    TokenKind.Rxn: NotImplementedFunc,
    TokenKind.BondDescriptor: map_bond_descriptor,
    TokenKind.StochasticSeperator: map_bond_seperator,
    TokenKind.StochasticStart: map_stochastic_object_start,
    TokenKind.StochasticEnd: map_stochastic_object_end,
    TokenKind.ImplictEndGroup: map_bond_descriptor,
    TokenKind.BondDescriptorLadder: NotImplementedFunc
}


def tokens_to_objects(constructor: BigSMILESStringConstructor, tokens: list[Token]) -> TokenKind | None:
    """

    Main loop for converting tokens into BigSMILES objects.

    Parameters
    ----------
    constructor: BigSMILESStringConstructor

    tokens:

    Returns
    -------

    """
    while tokens:
        token = tokens.pop(0)

        func = map_tokens[token.kind]
        result = func(constructor, tokens, token)
        if isinstance(result, TokenKind):
            return result

    constructor.run_validation()


def parse_bigsmiles_str(input_text: str, bigsmiles):
    """
    Main function that turns BigSMILES string into a BigSMILES object.
    Constructs BigSMILES tree in the provided object.
    """
    pre_validation(input_text)
    tokens = tokenize(input_text)

    constructor = BigSMILESStringConstructor(bigsmiles)
    tokens_to_objects(constructor, tokens)
    constructor.run_validation()
    constructor.exit()
