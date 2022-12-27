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
