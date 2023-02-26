import logging

from bigsmiles.errors import BigSMILESError
from bigsmiles.validation_string import run_string_validation
from bigsmiles.tokenizer import Token, TokenKind, tokenize
from bigsmiles.validation_tokens import run_token_validation
import bigsmiles.constructor_str as constructor
from bigsmiles.bigsmiles import BigSMILES, Branch, StochasticFragment, Atom, Bond, has_node_attr


def map_atom(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        return constructor.add_atom(parent, token.value)

    if isinstance(parent, StochasticFragment):
        return constructor.add_atom(parent, token.value)

    return constructor.add_bond_atom_pair(parent, "", token.value)


def map_bond(parent: has_node_attr, tokens: list[Token], token: Token):
    try:
        next_token = tokens.pop(0)
    except IndexError:
        raise BigSMILESError("Bond can't be at the end of a BigSMILES string.")

    if next_token.kind is TokenKind.Atom or next_token.kind is TokenKind.AtomExtend:
        return constructor.add_bond_atom_pair(parent, token.value, next_token.value)

    if next_token.kind is TokenKind.BondDescriptor:
        return constructor.add_bonding_descriptor_to_end_str(parent, token.value, next_token.value)

    if next_token.kind is TokenKind.StochasticStart:
        return map_stochastic_object_start(parent, tokens, next_token, token)

    raise BigSMILESError(f"Bond can't be followed by: {next_token.kind}")


def map_bond_descriptor(parent: has_node_attr, tokens: list[Token], token: Token):
    if tokens[0].kind is TokenKind.StochasticEnd:
        tokens.pop(0)
        parent = constructor.close_stochastic_fragment_str(parent)
        if token.kind is TokenKind.ImplictEndGroup:
            bond_symbol = None
        elif len(tokens) == 0 or tokens[0].kind is TokenKind.Atom or tokens[0].kind is TokenKind.AtomExtend:
            bond_symbol = ""
        elif tokens[0].kind is TokenKind.Bond:
            bond_symbol = tokens[0].value
        else:
            raise BigSMILESError("Bonding descriptors must be followed by an atom or bond.")

        return constructor.close_stochastic_object_str(parent, token.value, bond_symbol)

    if isinstance(parent, StochasticFragment):
        if not parent.nodes:
            # first symbol in stochastic fragment
            if tokens[0].kind == TokenKind.Atom:
                bond_symbol = ""
                atom_symbol = tokens.pop(0).value
            elif tokens[0].kind == TokenKind.Bond:
                bond_symbol = tokens[0].value
                tokens.pop(0)
                if tokens[0].kind == TokenKind.Atom:
                    atom_symbol = tokens.pop(0).value
                else:
                    raise BigSMILESError("Invalid SMILES. Bond descriptor followed by bond; "
                                         "must also be followed by atom.")
            else:
                raise BigSMILESError("Bonding descriptors must be followed by an atom or bond.")

            return constructor.add_bonding_descriptor_atom_pair_str(parent, atom_symbol, token.value, bond_symbol)

    elif isinstance(parent, Branch):
        # check to make sure the branch closes immediately
        if tokens[0].kind != TokenKind.BranchEnd:
            raise BigSMILESError("If Bonding Descriptors is the first Branch symbol, it can only be followed by "
                                 "Branch End.")

    return constructor.add_bonding_descriptor_to_end_str(parent, token.value, "")


def map_branch_start(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        raise BigSMILESError("BigSMILES can't start with a branch symbol.")
    if isinstance(parent, Branch) and len(parent.nodes) == 0:
        raise BigSMILESError("Branch can't be the first object within a branch.")
    return constructor.open_branch(parent)


def map_branch_end(parent: has_node_attr, tokens: list[Token], token: Token):
    return constructor.close_branch(parent)


def map_ring(parent: has_node_attr, tokens: list[Token], token: Token):
    return constructor.add_ring(parent, int(token.value.replace('%', '')))


def map_stochastic_object_start(parent: has_node_attr, tokens: list[Token], token: Token,
                                bond_token: Token = None):
    try:
        next_token = tokens.pop(0)
    except IndexError:
        raise BigSMILESError("Stochastic objects must begin with a bond descriptor (or implicit bonding descriptor).")

    if next_token.kind is TokenKind.ImplictEndGroup:
        return constructor.open_stochastic_object_str(parent, next_token.value, None)

    if next_token.kind is not TokenKind.BondDescriptor:
        raise BigSMILESError(f"Stochastic object starts must be followed an explict or implicit end group.")

    if isinstance(parent, BigSMILES) and not parent:
        return constructor.open_stochastic_object_str(parent, next_token.value, "")

    if bond_token is None:
        bond_token = ""
    return constructor.open_stochastic_object_str(parent, next_token.value, bond_token)


def map_stochastic_object_end(parent: has_node_attr, tokens: list[Token], token: Token) -> TokenKind:
    raise BigSMILESError("Stochastic objects should end with bonding descriptor (or implicit bonding description)")


def map_bond_seperator(parent: has_node_attr, tokens: list[Token], token: Token) -> has_node_attr:
    return constructor.close_open_stochastic_fragment_str(parent)


def NotImplementedFunc(*args, **kwargs):
    raise NotImplementedError()


def SkipSymbol(parent: has_node_attr, tokens: list[Token], token: Token) -> has_node_attr:
    logging.warning(f"Symbol skipped: {token.value}")
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


def tokens_to_bigsmiles(parent: has_node_attr, tokens: list[Token]):
    """

    Main loop for converting tokens into BigSMILES objects.

    Parameters
    ----------
    parent: has_node_attr

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
            raise BigSMILESError(f"Issue with token '{token}'. (token: {num_tokens-len(tokens)-1})", e) from e


def parse_bigsmiles_str(input_text: str, bigsmiles: BigSMILES):
    """
    Main function that turns BigSMILES string into a BigSMILES object.
    Constructs BigSMILES tree in the provided object.
    """
    input_text = run_string_validation(input_text)
    tokens = tokenize(input_text)
    tokens = run_token_validation(tokens)

    try:
        tokens_to_bigsmiles(bigsmiles, tokens)
        constructor.exit_construction(bigsmiles)
    except BigSMILESError as e:
        raise BigSMILESError(f"Parsing failed on '{input_text}'.", e) from e
