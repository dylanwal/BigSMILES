import warnings


def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'

warnings.formatwarning = custom_formatwarning

from bigsmiles.errors import BigSMILESError
from bigsmiles.validation import pre_validation
from bigsmiles.tokenizer import Token, TokenKind, tokenize
import bigsmiles.constructor_str as constructor
from bigsmiles.bigsmiles import BigSMILES, Branch, BondDescriptorAtom, Bond, Atom, StochasticObject, StochasticFragment


def map_atom(parent: constructor.ParentType, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        return constructor.add_atom(parent, token.value)

    if isinstance(parent, StochasticFragment):
        return constructor.add_atom(parent, token.value)

    return constructor.add_bond_atom_pair(parent, "", token.value)


def map_bond(parent: constructor.ParentType, tokens: list[Token], token: Token):
    try:
        next_token = tokens.pop(0)
    except IndexError:
        raise BigSMILESError("Bond can't be at the end of a BigSMILES string.")

    if next_token.kind in (TokenKind.Atom, TokenKind.AtomExtend):
        return constructor.add_bond_atom_pair(parent, token.value, next_token.value)

    if next_token.kind is not TokenKind.BondDescriptor:
        return constructor.add_bond_bonding_descriptor_pair_str(parent, token.value, next_token.value)

    return map_stochastic_object_start(parent, tokens, next_token, token)


def map_bond_descriptor(parent: constructor.ParentType, tokens: list[Token], token: Token):
    if tokens[0].kind is TokenKind.StochasticEnd:
        parent = constructor.close_stochastic_fragment_str(parent)
        parent = constructor.close_stochastic_object_str(parent, token.value)
        tokens.pop(0)
        return parent

    if isinstance(parent, StochasticFragment):
        # first StochasticFragment symbol
        return constructor.add_bonding_descriptor_str(parent, token.value)

    if isinstance(parent, BigSMILES) and not parent:
        # first Branch symbol
        # check to make sure the branch closes immediately
        if tokens[0].kind != TokenKind.BranchEnd:
            raise BigSMILESError("If Bonding Descriptors is the first Branch symbol, it can only be followed by "
                                 "Branch End.")

    return constructor.add_bond_bonding_descriptor_pair_str(parent, "", token.value)


def map_branch_start(parent: constructor.ParentType, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        raise BigSMILESError("BigSMILES can't start with a branch symbol.")
    if isinstance(parent, Branch) and len(parent.nodes) == 0:
        raise BigSMILESError("Branch can't be the first object within a branch.")
    return constructor.open_branch(parent)


def map_branch_end(parent: constructor.ParentType, tokens: list[Token], token: Token):
    return constructor.close_branch(parent)


def map_ring(parent: constructor.ParentType, tokens: list[Token], token: Token):
    return constructor.add_ring(parent, int(token.value.replace('%', '')))


def map_stochastic_object_start(parent: constructor.ParentType, tokens: list[Token], token: Token,
                                bond_token: Token = None):
    try:
        next_token = tokens.pop(0)
    except IndexError:
        raise BigSMILESError("Stochastic objects must begin with a bond descriptor (or implicit bonding descriptor).")

    if next_token.kind not in (TokenKind.ImplictEndGroup, TokenKind.BondDescriptor):
        raise BigSMILESError(f"Stochastic object starts must be followed an explict or implicit end group.")

    if isinstance(parent, BigSMILES) and not parent:
        return constructor.open_stochastic_object_fragment_str(parent, next_token.value)
    else:
        if bond_token is None:
            return constructor.open_stochastic_object_with_bond_str(parent, "", next_token.value)
        else:
            return constructor.open_stochastic_object_with_bond_str(parent, bond_token.value, next_token.value)


def map_stochastic_object_end(parent: constructor.ParentType, tokens: list[Token], token: Token) -> TokenKind:
    raise BigSMILESError("Stochastic objects should end with bonding descriptor (or implicit bonding description)")


def map_bond_seperator(parent: constructor.ParentType, tokens: list[Token], token: Token):
    return constructor.close_open_stochastic_fragment_str(parent)


def NotImplementedFunc(*args, **kwargs):
    raise NotImplementedError()


def SkipSymbol(parent: constructor.ParentType, tokens: list[Token], token: Token):
    warnings.warn(f"Symbol skipped: {token.value}")
    return parent


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


def tokens_to_bigsmiles(parent: constructor.ParentType, tokens: list[Token]):
    """

    Main loop for converting tokens into BigSMILES objects.

    Parameters
    ----------
    parent: constructor.ParentType

    tokens:

    Returns
    -------

    """
    num_tokens = len(tokens)
    while tokens:
        token = tokens.pop(0)
        func = map_tokens[token.kind]
        try:
            parent = func(parent, tokens, token)
        except BigSMILESError as e:
            raise BigSMILESError(f"Issue with token '{token}'. (index: {num_tokens-len(tokens)-1})", e) from e


def parse_bigsmiles_str(input_text: str, bigsmiles):
    """
    Main function that turns BigSMILES string into a BigSMILES object.
    Constructs BigSMILES tree in the provided object.
    """
    pre_validation(input_text)
    tokens = tokenize(input_text)

    try:
        tokens_to_bigsmiles(bigsmiles, tokens)
        constructor.run_validation(bigsmiles)
        constructor.exit_(bigsmiles)
    except BigSMILESError as e:
        raise BigSMILESError(f"Parsing failed on '{input_text}'.", e) from e
