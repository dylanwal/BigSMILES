
from bigsmiles.data_structures.bigsmiles import BigSMILES, Atom, Bond


class StochasticNode:
    def __init__(self, symbol: str, index: int, stoch_obj_id: int):
        self.symbol = symbol
        self.index = index
        self.stoch_obj_id = stoch_obj_id


class Edge:
    def __init__(self):
        pass


class DiGraph:
    def __init__(self, bigsmiles: BigSMILES):
        self.bigsmiles = bigsmiles

        self._edges = []
        self._nodes = []

    @property
    def edges(self):
        return self._edges

    @property
    def nodes(self):
        return self._nodes

    def has_node(self, node) -> bool:
        if node in self._nodes:
            return True

        return False

    def has_edge(self, edge) -> bool:
        if edge in self._edges:
            return True

        return False

    def add_node(self, node):
        self._nodes.append(node)

    def add_edge(self, edge: Bond):
        if edge.atom1 not in self._nodes:
            self.add_node(edge.atom1)
        if edge.atom2 not in self._nodes:
            self.add_node(edge.atom2)

        self._edges.append(edge)

    def remove_node(self, node):
        self._nodes.remove(node)
        node.delete()

