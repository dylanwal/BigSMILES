import logging

from bigsmiles.errors import BigSMILESError
from bigsmiles.validation.validation_string import run_string_validation
from bigsmiles.constructors.tokenizer import Token, TokenKind, tokenize
from bigsmiles.validation.validation_tokens import run_token_validation
import bigsmiles.constructors.constructor_str as constructor
from bigsmiles.data_structures.bigsmiles import BigSMILES, Branch, StochasticFragment, Atom, StochasticObject, has_node_attr


def map_atom(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        return constructor.add_atom_str(parent, token.value)

    if isinstance(parent, StochasticFragment) and len(parent.nodes) == 0:
        return constructor.add_atom_str(parent, token.value)

    if token.kind is TokenKind.Aromatic and constructor.get_prior(parent, (Atom, StochasticObject)).aromatic:
        return constructor.add_bond_atom_pair_str(parent, ":", token.value)

    return constructor.add_bond_atom_pair_str(parent, "", token.value)


def map_bond(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        raise BigSMILESError(f"Bond can't be the first symbol.")

    try:
        next_token = tokens.pop(0)
    except IndexError:
        raise BigSMILESError("Bond can't be at the end of a BigSMILES string.")

    if next_token.kind in (TokenKind.Atom, TokenKind.AtomExtend, TokenKind.Aromatic):
        return constructor.add_bond_atom_pair(parent, token.value, next_token.value)

    if next_token.kind is TokenKind.BondDescriptor:
        return constructor.add_bond_bonding_descriptor_pair_str(parent, token.value, next_token.value)

    if next_token.kind is TokenKind.StochasticStart:
        return map_stochastic_object_start(parent, tokens, next_token, token)

    if next_token.kind is TokenKind.Ring or next_token.kind is TokenKind.Ring2:
        return constructor.add_ring(parent, int(next_token.value.replace('%', '')), token.value)

    raise BigSMILESError(f"Bond can't be followed by: {next_token.kind}")


def map_bond_descriptor(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        raise BigSMILESError("Bonding Descriptors can't be the first symbol.")

    if tokens[0].kind is TokenKind.StochasticEnd:
        parent = constructor.close_stochastic_fragment_str(parent)
        parent = constructor.close_stochastic_object_str(parent, token.value)
        tokens.pop(0)
        return parent

    # if tokens[0].kind != TokenKind.BranchEnd:
    #     # check to make sure the branch closes immediately
    #     raise BigSMILESError("If Bonding Descriptors is the first Branch symbol, it can only be followed by "
    #                          "Branch End.")

    if isinstance(parent, StochasticFragment) and len(parent.nodes) == 0:
        # first StochasticFragment symbol
        return constructor.add_bonding_descriptor_str(parent, token.value)

    return constructor.add_bond_bonding_descriptor_pair_str(parent, "", token.value)


def map_branch_start(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        raise BigSMILESError("BigSMILES can't start with a branch symbol.")
    if isinstance(parent, Branch) and len(parent.nodes) == 0:
        raise BigSMILESError("Branch can't be the first object within a branch.")
    return constructor.open_branch(parent)


def map_branch_end(parent: has_node_attr, tokens: list[Token], token: Token):
    return constructor.close_branch(parent)


def map_ring(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        raise BigSMILESError(f"Ring can't be the first symbol.")

    return constructor.add_ring(parent, int(token.value.replace('%', '')))


def map_stochastic_object_start(parent: has_node_attr, tokens: list[Token], token: Token,
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


def map_stochastic_object_end(parent: has_node_attr, tokens: list[Token], token: Token) -> TokenKind:
    raise BigSMILESError("Stochastic objects should end with bonding descriptor (or implicit bonding description)")


def map_bond_seperator(parent: has_node_attr, tokens: list[Token], token: Token) -> has_node_attr:
    return constructor.close_open_stochastic_fragment_str(parent)


def map_reaction(parent: has_node_attr, tokens: list[Token], token: Token):
    raise BigSMILESError("Reaction Symbol detected. Use 'bigsmiles.Reaction()' to process this input.")


def NotImplementedFunc(*args, **kwargs):
    raise NotImplementedError()


def SkipSymbol(parent: has_node_attr, tokens: list[Token], token: Token) -> has_node_attr:
    logging.warning(f"Symbol skipped: {token.value}")
    return parent


# this dict maps a token kind to a function that will process or map the symbol into a BigSMILES
map_token_functions = {
    TokenKind.Bond: map_bond,
    TokenKind.Atom: map_atom,
    TokenKind.Aromatic: map_atom,
    TokenKind.AtomExtend: map_atom,
    TokenKind.BranchStart: map_branch_start,
    TokenKind.BranchEnd: map_branch_end,
    TokenKind.Ring: map_ring,
    TokenKind.Ring2: map_ring,
    TokenKind.BondEZ: SkipSymbol,
    TokenKind.Disconnected: NotImplementedFunc,
    TokenKind.Rxn: map_reaction,
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
    """
    num_tokens = len(tokens)
    while tokens:
        token = tokens.pop(0)
        func = map_token_functions[token.kind]
        try:
            parent = func(parent, tokens, token)
        except BigSMILESError as e:
            raise BigSMILESError(f"Issue with token '{token}'. (token: {num_tokens-len(tokens)-1})", e) from e


def parse_bigsmiles_str(input_text: str, bigsmiles: BigSMILES):
    """
    Main function that turns BigSMILES string tokens then into a BigSMILES object.
    """
    input_text = run_string_validation(input_text)
    tokens = tokenize(input_text)
    tokens = run_token_validation(tokens)

    try:
        tokens_to_bigsmiles(bigsmiles, tokens)
        constructor.exit_construction(bigsmiles)
    except BigSMILESError as e:
        raise BigSMILESError(f"Parsing failed on '{input_text}'.", e) from e
