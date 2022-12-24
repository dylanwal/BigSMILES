import types

try:
    import networkx as nx
except ImportError:
    raise ImportError("To use this feature install networkx. (pip install networkx)")

from bigsmiles.bigsmiles import BigSMILES, StochasticObject, StochasticFragment, Branch, Bond, BondDescriptor, BondDescriptorAtom, Atom
from bigsmiles.nx_graph.draw_nx_graph import draw


def create_nx_graph(bigsmiles: BigSMILES) -> nx.Graph:
    """ Starting point for creating nx.graph"""
    graph = nx.Graph(attr={'bigsmiles': bigsmiles})
    graph.draw = types.MethodType(draw, graph)

    for node in bigsmiles.nodes:
        add_obj(graph, node)

    return graph


remove_attr = {
    Atom: {"id_", "bonds"},
    Bond: {"atom1", "atom2"},
    BondDescriptor: {}
}


def get_obj_attr(obj: Atom | Bond | BondDescriptor) -> dict:
    """ Creates dict of attributes for object. """
    dict_ = {s: getattr(obj, s) for s in obj.__slots__ if hasattr(obj, s)}
    for attr in remove_attr[type(obj)]:
        if attr in dict_:
            del dict_[attr]

    return dict_


def add_atom(graph: nx.Graph, atom: Atom) -> str:
    label = f"{atom.symbol}{atom.id_}"
    if graph.has_node(label):
        return label

    graph.add_node(label, **get_obj_attr(atom))
    for bond in atom.bonds:
        add_bond(graph, bond)

    return label


def add_bond(graph: nx.Graph, bond: Bond):
    atom1_label = add_obj(graph, bond.atom1)
    atom2_label = add_obj(graph, bond.atom2)
    if graph.has_edge(atom1_label, atom2_label):
        return

    graph.add_edge(atom1_label, atom2_label, **get_obj_attr(bond))


def add_bonding_descriptor_atom(graph: nx.Graph, bond_descr: BondDescriptorAtom) -> str:
    return add_bonding_descriptor(graph, bond_descr.descriptor)


def add_bonding_descriptor(graph: nx.Graph, bond_descr: BondDescriptor) -> str:
    label = f"{bond_descr.symbol}{bond_descr.index_}" + "\n{}" + str(bond_descr.stochastic_object.id_)
    if graph.has_node(label):
        return label

    graph.add_node(label)

    return label


def add_branch(graph: nx.Graph, branch: Branch):
    for node in branch.nodes:
        add_obj(graph, node)


def add_stochastic_object(graph: nx.Graph, stoch_obj: StochasticObject) -> str:
    label = "{}" + str(stoch_obj.id_)
    if graph.has_node(label):
        return label

    graph.add_node(label)
    for bond_descr in stoch_obj.bonding_descriptors:
        bond_descr_label = add_bonding_descriptor(graph, bond_descr)
        graph.add_edge(label, bond_descr_label)

    for stoch_frag in stoch_obj.nodes:
        add_stochastic_fragment(graph, stoch_frag)

    return label


def add_stochastic_fragment(graph: nx.Graph, stoch_frag: StochasticFragment):
    for node in stoch_frag.nodes:
        add_obj(graph, node)


obj_func = {
    Atom: add_atom,
    StochasticObject: add_stochastic_object,
    StochasticFragment: add_stochastic_fragment,
    Branch: add_branch,
    BondDescriptorAtom: add_bonding_descriptor_atom,
    Bond: add_bond
}


def add_obj(graph, obj: Atom | StochasticObject) -> str:
    return obj_func[type(obj)](graph, obj)
