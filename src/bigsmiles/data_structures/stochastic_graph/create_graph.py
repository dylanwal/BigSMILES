import copy
import enum
import types

from bigsmiles.data_structures.stochastic_graph.graph import StochasticNode, Edge, DiGraph

from bigsmiles.data_structures.bigsmiles import BigSMILES, StochasticObject, StochasticFragment, Branch, Bond, \
    BondDescriptor, BondDescriptorAtom, Atom


def bigsmiles_to_stochastic_graph(bigsmiles: BigSMILES) -> DiGraph:
    """ Starting point for creating BigSMILES stochastic_graph. """
    graph = get_graph(bigsmiles)
    
    for i in range(10):
        stoch_objs = get_stochastic_objects(graph)
        if not stoch_objs:
            break

        for stoch_obj in stoch_objs:
            process_stochastic_objects(graph, stoch_obj)

    return graph


def get_graph(obj: BigSMILES | StochasticFragment) -> DiGraph:
    """ Creates a new Digraph and loops through nodes"""
    graph = DiGraph(obj)

    # add first path
    for node in obj.nodes:
        add_obj(graph, node)

    add_rings(graph, obj.rings)

    return graph


def get_stochastic_objects(graph: DiGraph) -> list[StochasticObject]:
    stochastic_objects = []
    for node in graph.nodes:
        if isinstance(node, StochasticObject):
            stochastic_objects.append(node)
    return stochastic_objects


def add_rings(graph: DiGraph, rings: list[Bond]):
    for ring in rings:
        add_bond(graph, ring)


def add_atom(graph: DiGraph, atom: Atom):
    graph.add_node(atom)


def add_bond(graph: DiGraph, bond: Bond):
    graph.add_edge(bond)


def add_branch(graph: DiGraph, branch: Branch):
    for node in branch.nodes:
        add_obj(graph, node)


def add_stochastic_object(graph: DiGraph, stoch_obj: StochasticObject):
    graph.add_node(stoch_obj)


def add_bonding_descriptor_atom(graph: DiGraph, bond_descr: BondDescriptorAtom):
    graph.add_node(bond_descr)


obj_func = {
    Atom: add_atom,
    StochasticObject: add_stochastic_object,
    Branch: add_branch,
    BondDescriptorAtom: add_bonding_descriptor_atom,
    Bond: add_bond
}


def add_obj(graph: DiGraph, obj: Atom | StochasticObject):
    obj_func[type(obj)](graph, obj)


#####################################################################################

class GraphBondDirection(enum.Enum):
    In = 0
    Out = 1


class GraphBond:
    def __init__(self, direction: GraphBondDirection, u: str, v: str, symbol: str):
        self.direction = direction
        self.u = u
        self.v = v
        self.symbol = symbol

    def __repr__(self):
        return f"{self.u} --> {self.v}"

    def reverse(self):
        self.u, self.v = self.v, self.u


def get_bonds_into_node(graph: DiGraph, node: str) -> list[GraphBond]:
    bonds_nodes = list(graph.in_edges(node))

    bonds = []
    for bond in bonds_nodes:
        bonds.append(GraphBond(GraphBondDirection.In, bond[0], bond[1], graph.get_edge_data(*bond)))

    return bonds


def get_bonds_out_of_node(graph: DiGraph, node: str) -> list[GraphBond]:
    bonds_nodes = list(graph.out_edges(node))

    bonds = []
    for bond in bonds_nodes:
        bonds.append(GraphBond(GraphBondDirection.Out, bond[0], bond[1], graph.get_edge_data(*bond)))

    return bonds


def get_all_bonds_of_node(graph: DiGraph, node: str) -> list[GraphBond]:
    bonds = []
    bonds += get_bonds_into_node(graph, node)
    bonds += get_bonds_out_of_node(graph, node)
    return bonds


def is_complement_bond_descr(bond_descr1: BondDescriptor, bond_descr2: BondDescriptor):
    if bond_descr1.index_ != bond_descr2.index_:
        return False

    if bond_descr1.descriptor == "$" and bond_descr2.descriptor == "$":
        return True
    if bond_descr1.descriptor == "<" and bond_descr2.descriptor == ">":
        return True
    if bond_descr1.descriptor == ">" and bond_descr2.descriptor == "<":
        return True

    return False


def process_stochastic_objects(graph: DiGraph, stoch_obj: StochasticObject):
    """ Main entry point for processing stochastic objects. """
    # grab end groups
    left_bond = GraphBond(stoch_obj.bond_left)
    right_bond = stoch_obj.bond_right
    graph.remove_node(stoch_obj)  # remove placeholder stochastic object

    # add
    bond_descr_dict = add_bond_descriptors(graph, stoch_obj)
    add_left_end_groups(graph, stoch_obj, left_bond, bond_descr_dict)
    add_right_end_groups(graph, stoch_obj, right_bond, bond_descr_dict)

    # add repeat units
    for stoch_frag in stoch_obj.nodes:
        add_stochastic_fragment(graph, stoch_frag)


def add_bond_descriptors(graph: DiGraph, stoch_obj: StochasticObject)-> dict:
    bond_descr_dict = {}
    for bond_descr in stoch_obj.bonding_descriptors:
        if any_unidirctional_stoch_frag(stoch_obj, bond_descr) or bond_descr.symbol == "$":
            node = StochasticNode(bond_descr.symbol, bond_descr.index_, bond_descr.parent.id_)
        else:
            node = StochasticNode("<>", bond_descr.index_, bond_descr.parent.id_)

        graph.add_node(node)
        bond_descr_dict[bond_descr] = node

    # attach bonding descriptor list to stochastic_graph for future use
    return bond_descr_dict


def add_left_end_groups(graph: DiGraph, stoch_obj: StochasticObject, left_bond: GraphBond, bond_descr_dict):
    for bd, bd_node in bond_descr_dict.items():
        if bd is stoch_obj.bd_left:
            edge = Edge(left_bond.u, bd_node)
            graph.add_edge(edge)


def add_right_end_groups(graph: DiGraph, stoch_obj: StochasticObject, right_bond: GraphBond, bond_descr_dict):
    for bd, bd_node in bond_descr_dict.items():
        if is_complement_bond_descr(bd, stoch_obj.bd_right):
            edge = Edge(bd_node, right_bond.v)
            graph.add_edge(edge)


def add_stochastic_fragment(graph: DiGraph, stoch_frag: StochasticFragment):
    stoch_frag_graph = get_graph(stoch_frag)

    counter = 0
    for node in stoch_frag_graph.nodes:
        if stoch_frag_graph.nodes[node]['node_type'] is BondDescriptorAtom:
            for bd, bd_node in graph.bond_descr.items():
                if bd is stoch_frag_graph.nodes[node]['descriptor']:
                    if bd_node[:2] == "<>" and bd.descriptor == "<":
                        continue

                    stoch_frag_graph2 = copy.copy(stoch_frag_graph)
                    map_labels = {node: node + "\n" + str(counter) for node in stoch_frag_graph2.nodes}
                    stoch_frag_graph2 = relabel_nodes(stoch_frag_graph2, map_labels)
                    add_stoch_fragment_single(graph, stoch_frag_graph2, node + "\n" + str(counter), bd_node)
                    counter += 1


def add_stoch_fragment_single(
        graph: DiGraph,
        stoch_frag_graph: DiGraph,
        node: str,
        bd_node: str
):
    start_bond = get_all_bonds_of_node(stoch_frag_graph, node)[0]
    stoch_frag_graph.remove_node(node)

    # attach to stochastic_graph
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
        # stochastic_graph.add_edge(start_bond.v, bd_node, **start_bond.attr)


def connect_ends(graph: DiGraph, node: str, prior_node: str = None):
    edges = get_bonds_into_node(graph, node)

    # flip chains to be outs
    edges_in = [edge for edge in edges if edge.direction is GraphBondDirection.In]
    flip_edges = get_subgraph_needing_direction_flip(graph, edges_in, prior_node)
    for edge in flip_edges:
        flip_edge(graph, edge.u, edge.v)

    edges = get_bonds_out_of_node(graph, node)
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


def get_subgraph_needing_direction_flip(
        graph: DiGraph,
        edges_in: list[GraphBond],
        prior_atom: str
) -> list[GraphBond]:
    edges = []
    for edge in edges_in:
        if edge.u is prior_atom:
            continue
        if edge.u not in (v for v in graph.bond_descr.values()):
            edges.append(edge)
            edges_ = get_all_bonds_of_node(graph, edge.u)
            edges_in_ = [edge_ for edge_ in edges_ if edge_.direction is GraphBondDirection.In]
            edges += get_subgraph_needing_direction_flip(graph, edges_in_, edge.u)

    return edges


def any_unidirctional_stoch_frag(stoch_obj: StochasticObject, bond_descr: BondDescriptor) -> bool:
    """ Return True if [<]R[<] or [>]R[>] is present. """
    for stoch_frag in stoch_obj.nodes:
        if bond_descr in stoch_frag.bonding_descriptors:
            count = 0
            for atom in stoch_frag.nodes:
                if isinstance(atom, BondDescriptorAtom) and atom.descriptor == bond_descr:
                    count += 1

            if count == 2:
                return True

    return False


def flip_edge(graph: DiGraph, node1: str, node2: str):
    if not graph.has_edge(node1, node2) and graph.has_edge(node2, node1):
        node1, node2 = node2, node1
    attrs = graph.get_edge_data(node1, node2)
    graph.remove_edge(node1, node2)
    graph.add_edge(node2, node1, **attrs)
