import logging

import bigsmiles.errors as errors
from bigsmiles.validation.validation_string import run_string_validation
from bigsmiles.constructors.tokenizer import Token, TokenKind, tokenize
from bigsmiles.validation.validation_tokens import run_token_validation
import bigsmiles.constructors.constructor_str as constructor
from bigsmiles.data_structures.bigsmiles import BigSMILES, Branch, StochasticFragment, Atom, StochasticObject, \
    has_node_attr, BondDescriptorAtom


def get_next_token(tokens: list[Token], text: str = "") -> Token:
    try:
        return tokens.pop(0)
    except IndexError:
        raise errors.BigSMILESError(f"{text} can't be at the end of a BigSMILES string.")


atom_tokens = {TokenKind.Atom, TokenKind.AtomExtend, TokenKind.Aromatic}
aromatic_acceptors = (Atom, StochasticObject, BondDescriptorAtom)
map_next_token_to_bond = {
    TokenKind.Aromatic: ":",
    TokenKind.Atom: "",
    TokenKind.AtomExtend: "",
    TokenKind.StochasticStart: "",
    TokenKind.Disconnected: ".",
    TokenKind.BondDescriptor: "",
}

## Functions for mapping tokens to BigSMILES ##
#######################################################################################################################
def map_atom(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, BigSMILES) and not parent:
        return constructor.add_atom_str(parent, token.value)

    if isinstance(parent, StochasticFragment) and len(parent.nodes) == 0:
        return constructor.add_atom_str(parent, token.value)

    if token.kind is TokenKind.Aromatic and constructor.get_prior(parent, aromatic_acceptors).aromatic:
        return constructor.add_bond_atom_pair_str(parent, ":", token.value)

    return constructor.add_bond_atom_pair_str(parent, "", token.value)


def map_bond(parent: has_node_attr, tokens: list[Token], token: Token):
    next_token = get_next_token(tokens, "Bond")

    if next_token.kind in atom_tokens:
        return constructor.add_bond_atom_pair_str(parent, token.value, next_token.value)

    if next_token.kind is TokenKind.BondDescriptor:
        return constructor.add_bond_bonding_descriptor_pair_str(parent, token.value, next_token.value)

    if next_token.kind is TokenKind.StochasticStart:
        return map_stochastic_object_start(parent, tokens, next_token, token)

    if next_token.kind is TokenKind.Ring or next_token.kind is TokenKind.Ring2:
        return constructor.add_ring(parent, int(next_token.value.replace('%', '')), token.value)

    raise errors.BigSMILESError(f"Bond can't be followed by: {next_token.kind}")


def map_bond_descriptor(parent: has_node_attr, tokens: list[Token], token: Token):
    next_token = get_next_token(tokens, "Bond Descriptor")
    if next_token.kind is TokenKind.StochasticEnd:
        return map_stochastic_object_end(parent, tokens, token)

    if isinstance(parent, Branch) and next_token.kind != TokenKind.BranchEnd:
        # check to make sure the branch closes immediately
        raise errors.BigSMILESError("If Bonding Descriptors is the first Branch symbol, "
                                    "it can only be followed by Branch End.")

    if isinstance(parent, StochasticFragment) and len(parent.nodes) == 0:
        # first StochasticFragment symbol
        # get bond
        if next_token.kind is TokenKind.Bond:
            bond = next_token.value
        elif next_token.kind in map_next_token_to_bond:
            bond = map_next_token_to_bond[next_token.kind]
        elif next_token.kind in (TokenKind.Ring, TokenKind.Ring2):
            bond = None
        else:
            raise errors.BigSMILESError(f"Invalid Symbol following stochastic object close.")
        tokens.insert(0, next_token)

        return constructor.add_bonding_descriptor_str(parent, token.value, bond)

    tokens.insert(0, next_token)
    if constructor.get_prior(parent, (Atom, StochasticObject, BondDescriptorAtom)).aromatic:
        bond = ":"
    else:
        bond = ""

    return constructor.add_bond_bonding_descriptor_pair_str(parent, bond, token.value)


def map_branch_start(parent: has_node_attr, tokens: list[Token], token: Token):
    if isinstance(parent, Branch) and len(parent.nodes) == 0:
        raise errors.BigSMILESError("Branch can't be the first object within a branch.")
    return constructor.open_branch(parent)


def map_branch_end(parent: has_node_attr, tokens: list[Token], token: Token):
    return constructor.close_branch(parent)


def map_ring(parent: has_node_attr, tokens: list[Token], token: Token):
    if constructor.get_prior(parent, (Atom, StochasticObject)).aromatic:
        bond = ":"
    else:
        bond = ""

    return constructor.add_ring(parent, int(token.value.replace('%', '')), bond)


def map_stochastic_object_start(parent: has_node_attr, tokens: list[Token], token: Token, bond_token: Token = None):
    next_token = get_next_token(tokens, "Stochastic objects start symbol")

    if next_token.kind is TokenKind.ImplictEndGroup:
        return constructor.open_stochastic_object_fragment_str(parent, next_token.value)

    if next_token.kind is not TokenKind.BondDescriptor:
        raise errors.BigSMILESError(f"Stochastic object starts must be followed an explict or implicit end group.")

    if isinstance(parent, BigSMILES) and not parent:
        parent = constructor.add_atom_str(parent, "[H]")

    if constructor.get_prior(parent, (Atom, StochasticObject)).aromatic:
        bond = ":"
    elif bond_token is None:
        bond = ""
    else:
        bond = bond_token.value

    return constructor.open_stochastic_object_with_bond_str(parent, bond, next_token.value)



def map_stochastic_object_end(parent: has_node_attr, tokens: list[Token], token: Token) -> TokenKind:
    parent = constructor.close_stochastic_fragment_str(parent)

    # get bond
    if len(tokens) < 1:
        bond = None
    elif tokens[0].kind is TokenKind.Bond:
        bond = tokens[0].value
    elif tokens[0].kind in map_next_token_to_bond:
        bond = map_next_token_to_bond[tokens[0].kind]
    elif tokens[0].kind in (TokenKind.Ring, TokenKind.Ring2):
        bond = None
    else:
        raise errors.BigSMILESError(f"Invalid Symbol following stochastic object close.")

    return constructor.close_stochastic_object_str(parent, token.value, bond)


def map_bond_seperator(parent: has_node_attr, tokens: list[Token], token: Token) -> has_node_attr:
    return constructor.close_open_stochastic_fragment_str(parent)


def map_reaction(parent: has_node_attr, tokens: list[Token], token: Token):
    raise errors.BigSMILESError("Reaction Symbol detected. Use 'bigsmiles.Reaction()' to process this input.")


def map_disconnect(parent: has_node_attr, tokens: list[Token], token: Token):
    next_token = get_next_token(tokens, "Disconnect")
    if next_token.kind in atom_tokens:
        return constructor.add_bond_atom_pair_str(parent, ".", next_token.value)

    if next_token.kind is TokenKind.StochasticStart:
        return constructor.open_stochastic_object_with_bond_str(parent, ".", next_token.value)

    if next_token.kind is TokenKind.BondDescriptor:
        return constructor.add_bond_bonding_descriptor_pair_str(parent, ".", next_token.value)

    raise errors.BigSMILESError(f"Disconnect can't be followed by '{next_token.kind.name}'. ")


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
    TokenKind.Disconnected: map_disconnect,
    TokenKind.Rxn: map_reaction,
    TokenKind.BondDescriptor: map_bond_descriptor,
    TokenKind.StochasticSeperator: map_bond_seperator,
    TokenKind.StochasticStart: map_stochastic_object_start,
    TokenKind.StochasticEnd: map_stochastic_object_end,
    TokenKind.ImplictEndGroup: map_bond_descriptor,
    TokenKind.BondDescriptorLadder: NotImplementedFunc
}

valid_first_symbols = {TokenKind.Atom, TokenKind.AtomExtend, TokenKind.Aromatic, TokenKind.StochasticStart}


def tokens_to_bigsmiles(parent: has_node_attr, tokens: list[Token]):
    """
    Main loop for converting tokens into BigSMILES objects.
    """
    num_tokens = len(tokens)

    # guard statements
    if num_tokens < 1:
        raise errors.BigSMILESError(f"No BigSMILES symbols detected.")
    if tokens[0].kind not in valid_first_symbols:
        raise errors.BigSMILESError(f"BigSMILES can't start with a '{tokens[0].kind.name}' symbol.")

    # main loop
    while tokens:
        token = tokens.pop(0)
        func = map_token_functions[token.kind]
        try:
            parent = func(parent, tokens, token)
        except errors.BigSMILESError as e:
            raise errors.BigSMILESError(f"Issue with token '{token}'. (token: {num_tokens - len(tokens) - 1})",
                                        e) from e


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
    except errors.BigSMILESError as e:
        raise errors.BigSMILESError(f"Parsing failed on '{input_text}'.", e) from e
