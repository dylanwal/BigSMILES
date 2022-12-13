

from bigsmiles.tokenizer import Token, TokenKind
from bigsmiles.bigsmiles import Atom, Bond, BondDescriptor, Branch, StochasticFragment, StochasticObject, \
    BigSMILES, BigSMILESError


class BigSMILESConstructor:
    pass



def get_prior_atom_for_bond(parent, flag: bool = False) -> Atom | None:
    if parent.nodes:
        if isinstance(parent.nodes[-1], (Atom, BondDescriptor, StochasticObject)):
            return parent.nodes[-1]
        if isinstance(parent.nodes[-1], Branch):
            if isinstance(parent.nodes[-2], (Atom, BondDescriptor, StochasticObject)):
                return parent.nodes[-2]
            if isinstance(parent.nodes[-2], Branch):
                if isinstance(parent.nodes[-3], (Atom, BondDescriptor, StochasticObject)):
                    return parent.nodes[-3]
    elif hasattr(parent, "parent"):
        if flag:
            raise BigSMILESError("Ring inside ring is not allowed without adding an atom.")
        return get_prior_atom_for_bond(parent.parent, flag=True)

    return None


def map_bond(parent: BigSMILES | Branch, tokens: list[Token], token: Token) -> Atom:
    try:
        next_token = tokens.pop(0)
    except IndexError:
        raise BigSMILESError("Bond can't be at the end of a BigSMILES string.")

    if next_token.kind != TokenKind.Atom and next_token.kind != TokenKind.AtomExtend:
        raise BigSMILESError(f"Bonds must be followed by an Atom. (prior atom: {parent.nodes[-1]} | bond to"
                             f"ken: {token.value} | next token: {next_token.value} ({next_token.kind}))")

    next_atom = Atom(symbol=next_token.value, id_=parent._get_atom_id())
    prior_atom = get_prior_atom_for_bond(parent)
    if prior_atom is None:
        raise BigSMILESError(f"No prior atom found for bond. (prior nodes:{parent.nodes})")

    bond = Bond(token.value, atom1=prior_atom, atom2=next_atom, id_=parent._get_bond_id())
    parent.nodes.append(bond)
    parent.nodes.append(next_atom)

    return next_atom


def map_atom(parent: BigSMILES | Branch, tokens: list[Token], token: Token) -> Atom:
    if isinstance(parent, BigSMILES) and len(parent.nodes) == 0:
        # first BigSMILES symbol
        atom = Atom(symbol=token.value, id_=parent._get_atom_id())
        parent.nodes.append(atom)
        return atom

    if isinstance(parent, Branch) and len(parent.nodes) == 0:
        prior_atom = get_prior_atom_for_bond(parent.parent)  # first Branch symbol
    else:
        prior_atom = get_prior_atom_for_bond(parent)

    atom = Atom(symbol=token.value, id_=parent._get_atom_id())
    bond = Bond('', atom1=prior_atom, atom2=atom, id_=parent._get_bond_id())  # add single bond
    parent.nodes.append(bond)
    parent.nodes.append(atom)

    return atom


def map_branch(parent: BigSMILES | Branch, tokens: list[Token], token: Token) -> Branch:
    if isinstance(parent, BigSMILES) and len(parent.nodes) == 0:
        raise BigSMILESError("BigSMILES string or branch can't start with branch")

    branch = Branch(parent, id_=parent._get_branch_id())
    result = tokens_to_objects(branch, tokens)
    if result is not TokenKind.BranchEnd:
        raise BigSMILESError(f"No branch closure. Branch: {branch}")
    parent.nodes.append(branch)

    return branch


def map_branch_end(parent: BigSMILES | Branch, tokens: list[Token], token: Token) -> TokenKind:
    return TokenKind.BranchEnd


def map_ring(parent: BigSMILES | Branch, tokens: list[Token], token: Token) -> Bond:
    if not isinstance(parent.nodes[-1], Atom):
        raise BigSMILESError(f"Ring number must follow atoms. (prior token: {parent.nodes[-1]}, "
                             f"current token: {token.value})")

    ring_id = int(token.value)
    # check if ring already exists
    ring = parent._get_ring(ring_id)
    if ring is None or ring.atom2 is not None:
        ring = Bond("", atom1=parent.nodes[-1], id_=parent._get_bond_id(), ring=ring_id)
        parent._add_ring(ring)
        return ring

    ring.atom2 = parent.nodes[-1]

    return ring


def map_bond_descriptor(parent: StochasticObject | StochasticFragment | Branch, tokens: list[Token], token: Token) \
        -> BondDescriptor | StochasticFragment | TokenKind:
    if not parent.in_stochastic_object:
        raise BigSMILESError("Bonding descriptor must be in stochastic object")

    if tokens[0].kind == TokenKind.StochasticEnd:
        if not isinstance(parent, StochasticObject):
            # this BD is the end of stochastic object
            tokens.insert(0, token)
            return TokenKind.StochasticEnd

        bond_descriptor = BondDescriptor(symbol=token.value, id_=parent._get_bond_descriptor_id())
        parent.end_group_right = bond_descriptor
        return bond_descriptor

    if isinstance(parent, StochasticObject):
        if len(parent.nodes) == 0 and parent.end_group_left is None:
            # first StochasticObject symbol
            bond_descriptor = BondDescriptor(symbol=token.value, id_=parent._get_bond_descriptor_id())
            parent.end_group_left = bond_descriptor
            return bond_descriptor
        else:
            stochastic_fragment = StochasticFragment(parent=parent, id_=parent._get_stochastic_fragment_id())
            bond_descriptor = BondDescriptor(symbol=token.value, id_=parent._get_bond_descriptor_id())
            stochastic_fragment.nodes.append(bond_descriptor)
            tokens_to_objects(stochastic_fragment, tokens)
            parent.nodes.append(stochastic_fragment)
            return stochastic_fragment

    if isinstance(parent, StochasticFragment) and len(parent.nodes) == 0:
        # first StochasticFragment symbol
        bond_descriptor = BondDescriptor(symbol=token.value, id_=parent._get_bond_descriptor_id())
        parent.nodes.append(bond_descriptor)
        return bond_descriptor

    if isinstance(parent, Branch) and len(parent.nodes) == 0:
        # first Branch symbol
        # check to make sure the branch closes immediately
        if tokens[0].kind != TokenKind.BranchEnd:
            raise BigSMILESError("If Bonding Descriptors is the first Branch symbol, it can only be followed by "
                                 "Branch End.")
        prior_atom = get_prior_atom_for_bond(parent.parent)
    else:
        prior_atom = get_prior_atom_for_bond(parent)

    bond_descriptor = BondDescriptor(symbol=token.value, id_=parent._get_bond_descriptor_id())
    bond = Bond('', atom1=prior_atom, atom2=bond_descriptor, id_=parent._get_bond_id())  # add single bond
    parent.nodes.append(bond)
    parent.nodes.append(bond_descriptor)

    return bond_descriptor


def map_stochastic_object(parent: BigSMILES | StochasticFragment | Branch, tokens: list[Token], token: Token) \
        -> StochasticObject:

    # check for prior bond
    if parent.nodes and not isinstance(parent.nodes[-1], Bond):
        prior_atom = get_prior_atom_for_bond(parent)
        stochastic_object = StochasticObject(parent, id_=parent._get_branch_id())
        bond = Bond('', atom1=prior_atom, atom2=stochastic_object, id_=parent._get_bond_id())  # add single bond
        parent.nodes.append(bond)

    else:
        stochastic_object = StochasticObject(parent, id_=parent._get_branch_id())

    result = tokens_to_objects(stochastic_object, tokens)
    if result is not TokenKind.StochasticEnd:
        raise BigSMILESError(f"No stochastic object closure. Stochastic_object: {stochastic_object}")
    parent.nodes.append(stochastic_object)

    return stochastic_object


def map_stochastic_object_end(parent: BigSMILES | Branch, tokens: list[Token], token: Token) -> TokenKind:
    return TokenKind.StochasticEnd


def map_bond_seperator(parent: StochasticFragment, tokens: list[Token], token: Token) -> TokenKind:
    if not isinstance(parent, StochasticFragment):
        raise BigSMILESError("Stochastic seperator can only follow a stochastic fragments.")

    return TokenKind.StochasticSeperator


def NotImplementedFunc(*args, **kwargs):
    raise NotImplementedError


map_tokens = {
    TokenKind.Bond: map_bond,
    TokenKind.Atom: map_atom,
    TokenKind.Aromatic: map_atom,
    TokenKind.AtomExtend: map_atom,
    TokenKind.BranchStart: map_branch,
    TokenKind.BranchEnd: map_branch_end,
    TokenKind.Ring: map_ring,
    TokenKind.Ring2: map_ring,
    TokenKind.BondEZ: NotImplementedFunc,
    TokenKind.Mix: NotImplementedFunc,
    TokenKind.Rxn: NotImplementedFunc,
    TokenKind.BondDescriptor: map_bond_descriptor,
    TokenKind.StochasticSeperator: map_bond_seperator,
    TokenKind.StochasticStart: map_stochastic_object,
    TokenKind.StochasticEnd: map_stochastic_object_end,
    TokenKind.ImplictEndGroup: NotImplementedFunc,
    TokenKind.BondDescriptorLadder: NotImplementedFunc
}


def tokens_to_objects(parent: BigSMILES | StochasticObject | StochasticFragment | Branch, tokens: list[Token]) -> \
        TokenKind | None:
    while tokens:
        token = tokens.pop(0)

        func = map_tokens[token.kind]
        result = func(parent, tokens, token)
        if isinstance(result, TokenKind):
            return result

    return None
