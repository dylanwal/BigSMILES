import copy
import enum
import types

try:
    import networkx as nx
except ImportError:
    raise ImportError("To use this feature install networkx. (pip install networkx)")

from bigsmiles.bigsmiles import BigSMILES, StochasticObject, StochasticFragment, Branch, Bond, BondDescriptor, \
    BondDescriptorAtom, Atom, BondDescriptorTypes
from bigsmiles.nx_graph.draw_nx_graph import draw


class BigSMILESGraphError(Exception):
    pass


class BigSMILESGraph(nx.DiGraph):


    def __init__(self, bigsmiles: BigSMILES, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._added = set()

    def draw(self):
        import bigsmiles.nx_graph.draw_nx_graph
        bigsmiles.nx_graph.draw_nx_graph.draw(self)


remove_attr = {
        Atom: {"id_", "bonds"},
        Bond: {"atom1", "atom2"},
        BondDescriptorAtom: {}
    }


def create_nx_graph(bigsmiles: BigSMILES) -> nx.Graph:
    """ Starting point for creating nx.graph"""
    graph = get_graph(bigsmiles)

    # process stochastic objects
    for i in range(10):
        stoch_objs = get_stochastic_objects(graph)
        if not stoch_objs:
            break

        for stoch_obj in stoch_objs:
            graph = process_stochastic_objects(graph, stoch_obj)

    return graph


def get_graph(obj: BigSMILES | StochasticFragment):
    graph = nx.DiGraph()
    graph.draw = types.MethodType(draw, graph)
    graph._added = set()

    # add first path
    for node in obj.nodes:
        add_obj(graph, node)

    add_rings(graph, obj.rings)

    return graph


def get_stochastic_objects(graph):
    stochastic_objects = []
    for node in graph.nodes:
        if "stoch_obj" in graph.nodes[node]:
            stochastic_objects.append(node)
    return stochastic_objects


def get_obj_attr(obj: Atom | Bond | BondDescriptorAtom) -> dict:
    """ Creates dict of attributes for object. """
    dict_ = {s: getattr(obj, s) for s in obj.__slots__ if hasattr(obj, s)}
    for attr in remove_attr[type(obj)]:
        if attr in dict_:
            del dict_[attr]

    return dict_


def add_rings(graph: nx.Graph, rings: list[Bond]):
    for ring in rings:
        add_bond(graph, ring)


def add_atom(graph: nx.Graph, atom: Atom) -> str:
    label = f"{atom.symbol}{atom.id_}"
    if graph.has_node(label):
        return label

    graph.add_node(label, node_type=Atom, **get_obj_attr(atom))
    return label


def add_bond(graph: nx.Graph, bond: Bond):
    atom1_label = add_obj(graph, bond.atom1)
    atom2_label = add_obj(graph, bond.atom2)

    if graph.has_edge(atom1_label, atom2_label):
        return

    graph.add_edge(atom1_label, atom2_label, **get_obj_attr(bond))


def add_branch(graph: nx.Graph, branch: Branch):
    for node in branch.nodes:
        add_obj(graph, node)


def add_stochastic_object(graph: nx.Graph, stoch_obj: StochasticObject) -> str:
    label = f"Q{stoch_obj.id_}"
    if graph.has_node(label):
        return label

    graph.add_node(label, node_type=StochasticObject, stoch_obj=stoch_obj)
    return label


def add_bonding_descriptor_atom(graph: nx.Graph, bond_descr: BondDescriptorAtom) -> str:
    label = f"{bond_descr.descriptor.symbol}{bond_descr.descriptor.index_}_{bond_descr.id_}" + \
            "\n{" + str(bond_descr.descriptor.stochastic_object.id_) + "}"
    if graph.has_node(label):
        return label

    graph.add_node(label, node_type=BondDescriptorAtom, **get_obj_attr(bond_descr))

    return label


obj_func = {
    Atom: add_atom,
    StochasticObject: add_stochastic_object,
    # StochasticFragment: add_stochastic_fragment,   # called from stochastic object
    Branch: add_branch,
    BondDescriptorAtom: add_bonding_descriptor_atom,
    Bond: add_bond
}


def add_obj(graph, obj: Atom | StochasticObject) -> str:
    return obj_func[type(obj)](graph, obj)


#####################################################################################

class GraphBondDirection(enum.Enum):
    In = 0
    Out = 1


class GraphBond:
    def __init__(self, direction: GraphBondDirection, u: str, v: str, attr: dict = None):
        self.direction = direction
        self.u = u
        self.v = v
        self.attr = attr

    def __repr__(self):
        return f"{self.u} --> {self.v}"

    def reverse(self):
        self.u, self.v = self.v, self.u


def get_bonds_of_node(graph: nx.DiGraph, node: str) -> list[GraphBond]:
    bonds = []
    bonds_in = list(graph.in_edges(node))
    for bond_ in bonds_in:
        bonds.append(GraphBond(GraphBondDirection.In, bond_[0], bond_[1], graph.get_edge_data(*bond_)))

    bonds_out = list(graph.out_edges(node))
    for bond_ in bonds_out:
        bonds.append(GraphBond(GraphBondDirection.Out, bond_[0], bond_[1], graph.get_edge_data(*bond_)))

    return bonds


def process_stochastic_objects(graph: nx.DiGraph, node: str):
    stoch_obj = graph.nodes[node]['stoch_obj']
    left_bond, right_bond = get_bonds_of_node(graph, node)

    graph.remove_node(node)
    sub_graphs = list(nx.weakly_connected_components(graph))

    left_graph = graph.subgraph(sub_graphs[0])
    right_graph = graph.subgraph(sub_graphs[1])

    subgraph = nx.DiGraph()
    subgraph.draw = types.MethodType(draw, subgraph)
    subgraph._added = set()

    # add bond descriptors
    bond_descr_dict = {}
    for bond_descr in stoch_obj.bonding_descriptors:
        if any_unidirctional_stoch_frag(stoch_obj, bond_descr) or bond_descr.symbol == "$":
            label = f"{bond_descr.symbol}{bond_descr.index_}" + "\n{}" + str(bond_descr.stochastic_object.id_)
            subgraph.add_node(label)
            bond_descr_dict[bond_descr] = label
            continue
        else:
            label = f"<>{bond_descr.index_}" + "\n{}" + str(bond_descr.stochastic_object.id_)
            subgraph.add_node(label)
            bond_descr_dict[bond_descr] = label
    subgraph.bond_descr = bond_descr_dict

    # Add left end group
    count = 0
    for bd, bd_node in subgraph.bond_descr.items():
        if bd is stoch_obj.end_group_left.descriptor:
            subgraph.add_edges_from(left_graph.edges(data=True))
            subgraph.add_edge(left_bond.u, bd_node, **left_bond.attr)
        if is_complement_bond_descr(bd, stoch_obj.end_group_left.descriptor):
            left_graph_2 = copy.deepcopy(left_graph)
            map_labels = {node: node + str(count) for node in left_graph_2.nodes}
            left_graph_2 = nx.relabel_nodes(left_graph_2, map_labels)
            subgraph.add_edges_from(left_graph_2.reverse().edges(data=True))
            subgraph.add_edge(bd_node, left_bond.u + str(count), **left_bond.attr)
            count += 1
            
    count = 0
    for bd, bd_node in subgraph.bond_descr.items():
        if bd is stoch_obj.end_group_right.descriptor:
            right_graph_2 = copy.deepcopy(right_graph)
            map_labels = {node: node + str(count) for node in right_graph_2.nodes}
            right_graph_2 = nx.relabel_nodes(right_graph_2, map_labels)
            subgraph.add_edges_from(right_graph_2.reverse().edges(data=True))
            subgraph.add_edge(right_bond.v + str(count), bd_node, **right_bond.attr)
            count += 1
        if is_complement_bond_descr(bd, stoch_obj.end_group_right.descriptor):
            subgraph.add_edges_from(right_graph.edges(data=True))
            subgraph.add_edge(bd_node, right_bond.v, **right_bond.attr)

    # add repeat units
    for stoch_frag in stoch_obj.nodes:
        add_stochastic_fragment(subgraph, stoch_frag)

    return subgraph


def is_complement_bond_descr(bond_descr1: BondDescriptor, bond_descr2: BondDescriptor):
    if bond_descr1.type_ is BondDescriptorTypes.Dollar and bond_descr2.type_ is BondDescriptorTypes.Dollar:
        return True
    if bond_descr1.type_ is BondDescriptorTypes.Left and bond_descr2.type_ is BondDescriptorTypes.Right:
        return True
    if bond_descr1.type_ is BondDescriptorTypes.Right and bond_descr2.type_ is BondDescriptorTypes.Left:
        return True

    return False


def add_stochastic_fragment(graph: nx.DiGraph, stoch_frag: StochasticFragment):
    stoch_frag_graph = get_graph(stoch_frag)

    counter = 0
    for node in stoch_frag_graph.nodes:
        if stoch_frag_graph.nodes[node]['node_type'] is BondDescriptorAtom:
            for bd, bd_node in graph.bond_descr.items():
                if bd is stoch_frag_graph.nodes[node]['descriptor']:
                    stoch_frag_graph2 = copy.copy(stoch_frag_graph)
                    map_labels = {node: node + "\n" + str(counter) for node in stoch_frag_graph2.nodes}
                    stoch_frag_graph2 = nx.relabel_nodes(stoch_frag_graph2, map_labels)
                    add_stoch_fragment_single(graph, stoch_frag_graph2, node + "\n" + str(counter), bd_node)
                    counter += 1


def add_stoch_fragment_single(
        graph: nx.DiGraph,
        stoch_frag_graph: nx.DiGraph,
        node: str,
        bd_node: str
):

    start_bond = get_bonds_of_node(stoch_frag_graph, node)[0]
    stoch_frag_graph.remove_node(node)

    # attach to graph
    if start_bond.direction == GraphBondDirection.Out:
        graph.add_nodes_from(stoch_frag_graph.nodes(data=True))
        graph.add_edges_from(stoch_frag_graph.edges(data=True))
        graph.add_edge(bd_node, start_bond.v, **start_bond.attr)
        connect_ends(graph, start_bond.v)

    else:
        # flip
        stoch_frag_graph = stoch_frag_graph.reverse()
        graph.add_nodes_from(stoch_frag_graph.nodes(data=True))
        graph.add_edges_from(stoch_frag_graph.edges(data=True))
        graph.add_edge(bd_node, start_bond.u, **start_bond.attr)
        connect_ends(graph, start_bond.u)
        # graph.add_edge(start_bond.v, bd_node, **start_bond.attr)


def connect_ends(graph: nx.DiGraph, node: str, prior_node: str = None):
    edges = get_bonds_of_node(graph, node)

    # flip chains to be outs
    edges_in = [edge for edge in edges if edge.direction is GraphBondDirection.In]
    flip_edges = flip_chain(graph, edges_in, prior_node)
    for edge in flip_edges:
        flip_edge(graph, edge.u, edge.v)

    edges = get_bonds_of_node(graph, node)
    edges_out = [edge for edge in edges if edge.direction is GraphBondDirection.Out]
    for edge in edges_out:
        if graph.nodes[edge.v]['node_type'] is BondDescriptorAtom:
            for bond_descr, bond_descr_node in graph.bond_descr.items():
                if is_complement_bond_descr(graph.nodes[edge.v]['descriptor'], bond_descr):
                    graph.add_edge(edge.u, bond_descr_node)
                    graph.remove_node(edge.v)
                    break
            else:
                raise BigSMILESGraphError('Code error.')
        else:
            connect_ends(graph, edge.v, node)


def flip_chain(graph: nx.DiGraph, edges_in: list[GraphBond], prior_atom: str) -> list[GraphBond]:
    edges = []
    for edge in edges_in:
        if edge.u is prior_atom:
            continue
        if edge.u not in (v for v in graph.bond_descr.values()):
            edges.append(edge)
            edges_ = get_bonds_of_node(graph, edge.u)
            edges_in_ = [edge_ for edge_ in edges_ if edge_.direction is GraphBondDirection.In]
            edges += flip_chain(graph, edges_in_, edge.u)

    return edges


def any_unidirctional_stoch_frag(stoch_obj: StochasticObject, bond_descr: BondDescriptor) -> bool:
    for stoch_frag in stoch_obj.nodes:
        if bond_descr in stoch_frag.bonding_descriptors:
            count = 0
            for atom in stoch_frag.nodes:
                if isinstance(atom, BondDescriptorAtom) and atom.descriptor == bond_descr:
                    count += 1

            if count == 2:
                return True

    return False


def flip_edge(graph: nx.DiGraph, node1: str, node2: str):
    if not graph.has_edge(node1, node2) and graph.has_edge(node2, node1):
        node1, node2 = node2, node1
    attrs = graph.get_edge_data(node1, node2)
    graph.remove_edge(node1, node2)
    graph.add_edge(node2, node1, **attrs)
#
#
# def get_all_simple_paths(subgraph: nx.DiGraph) -> list[list[str]]:
#     bond_descr_nodes = [node for node in subgraph.nodes if node[0] == "Q" or node[0] == "J"]
#
#     # if one node
#     if len(bond_descr_nodes) == 1:
#         cycles = get_cycles_with_bonding_descriptor(subgraph)
#         return [list(reversed(cycle)) for cycle in cycles]
#     else:
#         pass
#
#     return list(nx.all_simple_paths(subgraph, 'C2', 'C2'))
#
#
# def get_cycles_with_bonding_descriptor(subgraph) -> list[list[str]]:
#     cycles = nx.cycle_basis(subgraph.to_undirected())
#     cycles_with_bd = []
#     for cycle in cycles:
#         for node in cycle:
#             if node[0] == 'Q' or node[0] == 'J':
#                 cycles_with_bd.append(cycle)
#                 break
#
#     return cycles_with_bd
