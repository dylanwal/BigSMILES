import networkx as nx
import re
import copy
import rdkit
from rdkit import Chem
# import pydot
from networkx.algorithms.components.connected import connected_components


def getComp(d):
    if "<" in d:
        return d.replace("<", ">")
    elif ">" in d:
        return d.replace(">", "<")
    else:
        return d


def get_objects(ends):
    i = 0
    o = []
    start = []
    while i < len(ends):
        if ends[i] == "{":
            start.append(i)
            object = "{"
            count = 1
            i += 1
            while count != 0:
                if ends[i] == "{":
                    count += 1
                if ends[i] == "}":
                    count -= 1
                object += ends[i]
                i += 1
            o.append(object)
        else:
            i += 1
    return [o, start]


def get_repeats(object):
    object = object[object.find("]") + 1:]
    object = object[:object.rfind("[")]

    i = 0
    o = []
    count = 0
    repeat = ""
    object += ","
    while i < len(object):
        if object[i] == "{":
            count += 1
        if object[i] == "}":
            count -= 1
        if object[i] == "," and count == 0:
            o.append(repeat)
            repeat = ""
        else:
            repeat += object[i]
        i += 1

    repeats = [[], []]
    index = 0
    for i in range(len(o)):
        if ";" in o[i]:
            o[i].split(";")
            repeats[index].append(o[i][:o[i].index(";")])
            index += 1
            repeats[index].append(o[i][o[i].index(";") + 1:])
        else:
            repeats[index].append(o[i])
    return repeats[0], repeats[1]


def replace_objects_Bk(bigsmiles):
    smiles = ""
    counter = 0
    objects = get_objects(bigsmiles)[0]
    indices = get_objects(bigsmiles)[1]
    for i in range(len(objects)):
        smiles += bigsmiles[counter:indices[i]]
        smiles += "[Bk]"
        counter += (indices[i] - counter) + len(objects[i])
    smiles += bigsmiles[counter:]
    return smiles, objects


def add_Es(smiles, object_list):
    Bk_locations = [(d.start(0), d.end(0)) for d in re.finditer(r"\[Bk\]", smiles)]
    for i in range(len(Bk_locations)):
        repeats, implicit_ends = get_repeats(object_list[i])
        left_terminal = object_list[i][1:object_list[i].find("]") + 1]
        right_terminal = object_list[i][object_list[i].rfind("["):-1]
        repeat_unit_list = left_terminal
        for r in repeats:
            repeat_unit_list += r + ","
        repeat_unit_list += right_terminal

        if Bk_locations[i][0] == 0:
            descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(r"\[.\d+\]", repeat_unit_list)]
            a = descriptor_locations[-1][0]
            b = descriptor_locations[-1][1]
            if object_list[i].find("{[]") == 0:
                object_list[i] = "{" + getComp(repeat_unit_list[a:b]) + object_list[i][3:]
        if Bk_locations[i][1] == len(smiles) or smiles[Bk_locations[i][1]] == ")":
            descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(r"\[.\d+\]", repeat_unit_list)]
            a = descriptor_locations[0][0]
            b = descriptor_locations[0][1]
            if object_list[i][-3:] == "[]}":
                object_list[i] = object_list[i][0:-3] + getComp(repeat_unit_list[a:b]) + "}"

    while True:
        smiles_prev = smiles.replace("[Bk][Bk]", "[Bk][Es][Bk]")
        if smiles_prev == smiles:
            break
        smiles = smiles_prev
    while True:
        smiles_prev = smiles.replace("[Bk])", "[Bk][Es])")
        if smiles_prev == smiles:
            break
        smiles = smiles_prev
    while True:
        smiles_prev = smiles.replace("[Cf][Bk]", "[Cf][Es][Bk]")
        if smiles_prev == smiles:
            break
        smiles = smiles_prev
    while True:
        smiles_prev = smiles.replace("[Bk][Cf]", "[Bk][Es][Cf]")
        if smiles_prev == smiles:
            break
        smiles = smiles_prev

    if smiles.find("[Bk]") == 0:
        smiles = "[Es]" + smiles
    if smiles.rfind("[Bk]") == len(smiles) - 4 and smiles.rfind("[Bk]") != -1:
        smiles = smiles + "[Es]"

    return smiles, object_list


def RDKit_to_networkx_graph(mol):
    # https://github.com/maxhodak/keras-molecules/pull/32/files
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   symbol=atom.GetSymbol(),
                   formal_charge=atom.GetFormalCharge(),
                   is_aromatic=atom.GetIsAromatic(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   num_implicit_hs=atom.GetNumImplicitHs(),
                   active=False)
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
    return G


def orientation(networkx_graph, index):
    neighbors = []
    symbols = nx.get_node_attributes(networkx_graph, "symbol")
    i = -1
    for key in symbols:
        if symbols[key] == "Cf":
            i += 1
        if i == index:
            source = key
            break
    for t in nx.bfs_edges(networkx_graph, source=source):
        if symbols[t[1]] == "Bk":
            if t[1] > t[0]:
                neighbors.append([t[1], "1"])
            else:
                neighbors.append([t[1], "2"])
    neighbors = sorted(neighbors, key=lambda x: x[0])
    return [n[1] for n in neighbors]


def build_atomistic_directed(bigsmiles):
    bigsmiles = bigsmiles.replace("[<]", "[<1]").replace("[>]", "[>1]").replace("[$]", "[$1]")
    smiles, object_list = replace_objects_Bk(bigsmiles)
    smiles, object_list = add_Es(smiles, object_list)
    rdkit_graph = Chem.MolFromSmiles(smiles)
    networkx_graph = RDKit_to_networkx_graph(rdkit_graph)
    return build_level(networkx_graph, object_list, ["1"] * len(object_list))


def build_level(graph, object_list, object_orientation):
    symbols = nx.get_node_attributes(graph, "symbol")
    index_Bk = []
    for key in symbols:
        if symbols[key] == "Bk":
            neighbors = list(graph[key])
            if len(neighbors) == 2:
                index_Bk.append(key)

    nested_objects = []
    nested_orientation = []
    for o in range(len(object_list)):
        repeats, implicit_ends = get_repeats(object_list[o])
        left_terminal = object_list[o][1:object_list[o].find("]") + 1]
        right_terminal = object_list[o][object_list[o].rfind("["):-1]
        single_path = single_path_chemistries(repeats)

        neighbors = list(graph[index_Bk[o]])
        for n in neighbors:
            graph.remove_edge(index_Bk[o], n)

        def insert_terminals(graph, terminals, neighbors, bond):
            index = int(bond) - 1
            left_terminal = terminals[index][0]
            right_terminal = terminals[index][1]
            l_index = graph.number_of_nodes()
            graph.add_node(l_index, symbol=left_terminal, active=True)
            if left_terminal == right_terminal:
                r_index = l_index
            else:
                r_index = graph.number_of_nodes()
                graph.add_node(r_index, symbol=right_terminal, active=True)
            if bond == "1":
                graph.add_edge(neighbors[0], l_index, bond_type="1")
                graph.add_edge(neighbors[1], r_index, bond_type="2")
                terminals_attachment = [["2"], ["1"]]
            else:
                graph.add_edge(neighbors[0], l_index, bond_type="2")
                graph.add_edge(neighbors[1], r_index, bond_type="1")
                terminals_attachment = [["1"], ["2"]]
            return graph, terminals_attachment

        # add end group
        terminals = [[left_terminal, getComp(right_terminal)], [getComp(left_terminal), right_terminal]]
        if len(single_path["2"]) == 1:
            value = single_path["2"][0]
            if value == getComp(left_terminal):
                graph, terminals_attachment = insert_terminals(graph, terminals, neighbors, "1")
            elif value == getComp(right_terminal):
                graph, terminals_attachment = insert_terminals(graph, terminals, neighbors, "2")
        elif len(single_path["2"]) > 1:
            for value in single_path["2"]:
                if object_orientation[o] == "1":
                    if value == getComp(left_terminal):
                        graph, terminals_attachment = insert_terminals(graph, terminals, neighbors, "1")
                        single_path["2"] = [value]
                        break
                else:
                    if value == getComp(right_terminal):
                        graph, terminals_attachment = insert_terminals(graph, terminals, neighbors, "2")
                        single_path["2"] = [value]
                        break
        else:
            graph, terminals_attachment = insert_terminals(graph, terminals, neighbors, object_orientation[o])

        # for adding implicit endgorup
        if len(implicit_ends) > 0:
            symbols = nx.get_node_attributes(graph, "symbol")
            for i in range(2):
                n = list(graph[neighbors[i]])
                expl_a = symbols[neighbors[i]] != "Es"
                expl_b = symbols[neighbors[i]] == "Es" and len(list(graph[neighbors[i]])) >= 2
                if not (expl_a or expl_b):
                    terminals_attachment[i] = ["1", "2"]
            if left_terminal == getComp(right_terminal):
                a = list(set(terminals_attachment[0]) & set(terminals_attachment[1]))
                terminals_attachment = [a, a]


        # main loop for adding repeat units
        for smiles in repeats:
            smiles, nested_object_list = replace_objects_Bk(smiles)

            descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(r"\[.\d+\]", smiles)]
            descriptors = []
            for d in descriptor_locations:
                descriptors.append(smiles[d[0]:d[1]])
            for d in descriptors:
                smiles = smiles.replace(d, "[Cf]")

            duplication = len(descriptors)
            for i in range(len(descriptors)):
                if descriptors[i] in single_path["1"]:
                    duplication -= 1
                if getComp(descriptors[i]) in single_path["2"]:
                    duplication -= 1

            if "[Bk]" in smiles:
                smiles, nested_object_list = add_Es(smiles, nested_object_list)
                for d in range(duplication):
                    for n in nested_object_list:
                        nested_objects.append(n)
            rdkit_graph = Chem.MolFromSmiles(smiles)
            networkx_graph = RDKit_to_networkx_graph(rdkit_graph)

            for d in range(duplication):
                graph = nx.disjoint_union(graph, networkx_graph)

            symbols = nx.get_node_attributes(graph, "symbol")
            neighbors = []
            for key in symbols:
                if symbols[key] == "Cf":
                    n = list(graph[key])
                    if len(n) == 1:
                        neighbors.append(n[0])
            for key in symbols:
                if symbols[key] == "Cf":
                    n = list(graph[key])
                    if len(n) == 1:
                        graph.remove_edge(key, n[0])

            n = 0
            for i in range(len(descriptors)):
                if descriptors[i] in single_path["1"] or getComp(descriptors[i]) in single_path["2"]:
                    continue
                for o in range(len(descriptors)):
                    symbols = nx.get_node_attributes(graph, "symbol")
                    active = nx.get_node_attributes(graph, "active")
                    if i == o:
                        input = descriptors[i]
                        for key in symbols:
                            if symbols[key] == getComp(input) and active[key]:
                                graph.add_edge(key, neighbors[n], bond_type="2")
                                x = orientation(networkx_graph, i)
                                for j in x:
                                    nested_orientation.append(j)
                                break
                            elif key == len(symbols) - 1:
                                added = graph.number_of_nodes()
                                graph.add_node(added, symbol=getComp(input), active=True)
                                graph.add_edge(added, neighbors[n], bond_type="2")
                                x = orientation(networkx_graph, i)
                                for j in x:
                                    nested_orientation.append(j)
                                break
                    else:
                        output = descriptors[o]
                        for key in symbols:
                            if symbols[key] == output and active[key]:
                                graph.add_edge(key, neighbors[n], bond_type="1")
                                break
                            elif key == len(symbols) - 1:
                                added = graph.number_of_nodes()
                                graph.add_node(added, symbol=output, active=True)
                                graph.add_edge(added, neighbors[n], bond_type="1")
                                break
                    n += 1

        if len(implicit_ends) > 0:
            symbols = nx.get_node_attributes(graph, "symbol")
            active = nx.get_node_attributes(graph, "active")
            junctions = []
            for key in symbols:
                if active[key]:
                    descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(r"\[.\d+\]", symbols[key])]
                    if len(descriptor_locations) == 1:
                        junctions.append(key)

        for smiles in implicit_ends:
            smiles, nested_object_list = replace_objects_Bk(smiles)

            descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(r"\[.\d+\]", smiles)]
            descriptor = smiles[descriptor_locations[0][0]:descriptor_locations[0][1]]
            smiles = smiles.replace(descriptor, "[Cf]")
            if "[Bk]" in smiles:
                smiles, nested_object_list = add_Es(smiles, nested_object_list)
            rdkit_graph = Chem.MolFromSmiles(smiles)
            networkx_graph = RDKit_to_networkx_graph(rdkit_graph)

            def allowed_to_add(graph_descriptor, implicit_descriptor, bond_type):
                if bond_type == "1" and graph_descriptor == implicit_descriptor:
                    return True
                if bond_type == "2" and getComp(graph_descriptor) == implicit_descriptor:
                    return True
                return False

            def add_terminal(graph, smiles, graph_key, b_type):
                graph = nx.disjoint_union(graph, smiles)
                symbols = nx.get_node_attributes(graph, "symbol")
                for key in symbols:
                    if symbols[key] == "Cf":
                        n = list(graph[key])
                        if len(n) == 1:
                            neighbor = n[0]
                            graph.remove_edge(key, neighbor)
                            graph.add_edge(graph_key, neighbor, bond_type=b_type)
                return graph

            for key in junctions:
                if symbols[key] == left_terminal:
                    iteration = terminals_attachment[0]
                elif symbols[key] == right_terminal:
                    iteration = terminals_attachment[1]
                else:
                    iteration = ["1", "2"]
                for b_type in iteration:
                    if allowed_to_add(symbols[key], descriptor, b_type):
                        graph = add_terminal(graph, networkx_graph, key, b_type)
                        for n in nested_object_list:
                            nested_objects.append(n)
                            nested_orientation.append("1")

        nx.set_node_attributes(graph, False, "active")

    if len(nested_objects) > 0:
        graph = build_level(graph, nested_objects, nested_orientation)

    root = 0
    symbols = nx.get_node_attributes(graph, "symbol")
    for key in symbols:
        if key == root:
            continue
        x = list(nx.all_simple_paths(graph, root, key))
        if len(x) == 0:
            graph.remove_node(key)
        else:
            d = list(re.finditer(r"\[.\d+\]", symbols[key]))
            if len(d) > 0:
                neighbors = list(graph[key])
                if len(neighbors) == 2:
                    graph.remove_node(key)
                    graph.add_edge(*tuple(neighbors), bond_type=rdkit.Chem.rdchem.BondType.SINGLE)

    node_id = dict()
    for key in symbols:
        node_id[key] = str(key) + ": " + symbols[key]
    nx.set_node_attributes(graph, node_id, "node_id")

    return graph


def build_SMILES(graph):
    def extract_atoms(graph, extracted, start_atom, descriptors):
        extracted.add(start_atom)
        neighbors = list(graph[start_atom])
        next_atoms = []
        for n in neighbors:
            if n not in descriptors and n not in extracted:
                extracted.add(n)
                next_atoms.append(n)
        for n in next_atoms:
            extracted = extract_atoms(graph, extracted, n, descriptors)
        return extracted

    elements = nx.get_node_attributes(graph, "symbol")

    descriptors = []
    for key in elements:
        d = list(re.finditer(r"\[.\d+\]", elements[key]))
        if len(d) != 0:
            descriptors.append(key)

    edge_labels = nx.get_edge_attributes(graph, "bond_type")
    start = []
    for i in edge_labels:
        if edge_labels[i] == "1":
            d = list(re.finditer(r"\[.\d+\]", elements[i[0]]))
            if len(d) == 0:
                start.append(i[0])
            else:
                start.append(i[1])
    if len(start) == 0:
        start = [0]

    extracted = []
    for i in range(len(start)):
        extracted.append(extract_atoms(graph, set(), start[i], descriptors))
    ids = dict()
    counter = 1
    for i in extracted:
        for k in i:
            ids[k] = counter
        counter += 1
    nx.set_node_attributes(graph, ids, "ids")

    start = []
    for i in edge_labels:
        if edge_labels[i] == "2":
            d = list(re.finditer(r"\[.\d+\]", elements[i[0]]))
            if len(d) == 0:
                if i[0] not in ids:
                    start.append(i[0])
            else:
                if i[1] not in ids:
                    start.append(i[1])

    extracted = []
    for i in range(len(start)):
        extracted.append(extract_atoms(graph, set(), start[i], descriptors))
    for i in extracted:
        for k in i:
            ids[k] = counter
        counter += 1
    for i in descriptors:
        ids[i] = counter
        counter += 1
    nx.set_node_attributes(graph, ids, "ids")

    graph_directed = graph.copy()
    edge_labels = nx.get_edge_attributes(graph_directed, "bond_type")
    for d in descriptors:
        neighbors = graph_directed[d]
        same_frag = dict()
        for n in neighbors:
            same_frag[n] = ids[n]
        for i in same_frag:
            for j in same_frag:
                if i != j and same_frag[i] == same_frag[j]:
                    a = edge_labels[(min(d, i), max(d, i))]
                    b = edge_labels[(min(d, j), max(d, j))]
                    if [a, b] == ["1", "2"] or [b, a] == ["1", "2"]:
                        edge_labels[(min(d, i), max(d, i))] = "3"
                        edge_labels[(min(d, j), max(d, j))] = "3"
    nx.set_edge_attributes(graph_directed, edge_labels, "bond_type")

    edge_labels = nx.get_edge_attributes(graph_directed, "bond_type")
    vals = list(edge_labels.values())
    while vals.count("1") + vals.count("2") + vals.count("3") < len(vals):
        for i in edge_labels:
            if edge_labels[i] not in ["1", "2", "3"]:
                graph_directed = nx.contracted_edge(graph_directed, i, self_loops=False)
                break
        edge_labels = nx.get_edge_attributes(graph_directed, "bond_type")
        vals = list(edge_labels.values())

    graph_contracted = graph_directed
    graph_directed = nx.to_directed(graph_directed)
    graph_directed = nx.DiGraph(graph_directed)
    edge_labels = nx.get_edge_attributes(graph_directed, "bond_type")
    for i in edge_labels:
        a = edge_labels[i] == "1" and i[0] in descriptors
        b = edge_labels[i] == "2" and i[1] in descriptors
        if a or b:
            graph_directed.remove_edge(*i)
    edge_labels = nx.get_edge_attributes(graph_directed, "bond_type")
    for i in edge_labels:
        edge_labels[i] = ""
    nx.set_edge_attributes(graph_directed, edge_labels, "bond_type")

    return graph_contracted, graph_directed, descriptors


def single_path_chemistries(repeats):
    repeats = copy.deepcopy(repeats)
    forced = {"1": [], "2": []}

    # check for no complemntary bondinf descriptor
    for i in range(len(repeats)):
        descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(r"\[.\d+\]", repeats[i])]
        for d in descriptor_locations:
            desc1 = repeats[i][d[0]:d[1]]
            if "$" in desc1:
                continue
            no_compatible = True
            for j in range(len(repeats)):
                descriptor_locations = [(l.start(0), l.end(0)) for l in re.finditer(r"\[.\d+\]", repeats[j])]
                for l in descriptor_locations:
                    desc2 = repeats[j][l[0]:l[1]]
                    if desc1 == getComp(desc2):
                        no_compatible = False
            if no_compatible:
                forced["1"].append(desc1)

    descriptors = []
    for i in range(len(repeats)):
        for value in forced["1"]:
            repeats[i] = repeats[i].replace(value, "[]")
        descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(r"\[.\d+\]", repeats[i])]
        descriptors.append([])
        for d in descriptor_locations:
            desc1 = repeats[i][d[0]:d[1]]
            if "$" in desc1:
                return forced
            descriptors[-1].append(desc1)
    try:
        head = [descriptors[0][0], getComp(descriptors[0][0])]
        for i in range(len(head)):
            h = head[i]
            t = getComp(h)
            found_dendrimer = True
            for j in range(len(descriptors)):
                a = descriptors[j].count(h) + descriptors[j].count(t) == len(descriptors[j]) and descriptors[j].count(
                    h) == 1
                if not a:
                    found_dendrimer = False
                    break
            if found_dendrimer:
                forced["2"].append(h)
        return forced
    except:
        return forced


def visualize(graph):
    for g in graph:
        visual = nx.nx_agraph.write_dot(g, "visualize/dot.txt")
        text_file = open("visualize/dot.txt", "r")
        data = text_file.read()
        if "atomistic" in graph[g]:
            data = data.replace("node_id=", "shape=circle, label=")
        else:
            data = data.replace("ids=", "shape=circle, label=")
        data = data.replace("bond_type", "label")
        text_file.close()
        graphs = pydot.graph_from_dot_data(data)
        graphs[0].write_svg("visualize/" + graph[g] + ".svg")

# object_list = ["{[][>1]CCO[<1][>1]}","{[>1][<1]CCO[>1][<1]}","{[>1][>1]CCO[<1][]}","{[>1][>1]CCO[<1][]}"]
# smiles = "[Bk]CCC(CCC[Bk][Bk])[Bk]"
# smiles, object_list = add_Es(smiles, object_list)
# print(smiles)
# print(object_list)
# # print(smiles)
# # print(object_list)

# smiles = "[Cf]C[Bk]C[Bk]C[Bk]C[Cf]"
# RDKit_graph = Chem.MolFromSmiles(smiles)
# networkx_graph = RDKit_to_networkx_graph(RDKit_graph)
# x = orientation(networkx_graph, 1)
# print(x)

# query = "N{[<1][>1]CC(c1ccccc1)C\C=C/C[<1][>1]}N"
# query_atomistic_directed = build_atomistic_directed(query)
# query_contracted, query_SMILES_directed = build_SMILES(query_atomistic_directed)
# visualize({query_atomistic_directed: "Qatomistic_directed",
# query_SMILES_directed: "QSMILES_directed"})


def draw(self: nx.Graph):
    import matplotlib.pyplot as plt
    plt.close()
    plt.ylim([-1.1, 1.1])
    plt.xlim([-1.1, 1.1])
    plt.axis('off')
    nx.draw_networkx(self, pos=nx.kamada_kawai_layout(self),
                     node_size=100,
                     node_color=(0,0,0),
                     labels={node:self.nodes[node]['symbol']  for node in self.nodes},
                     width=5,
                     with_labels=True,
                     font_size=7,
                     font_color="w",
                     font_family="Arial Rounded MT Bold"
                     )
    plt.show()


def time_bigsmiles_parsing_graph(polymers: list[str], iter_: int) -> float:
    import datetime
    import time
    print(f"Starting graph time test: {datetime.datetime.now().time()}  (May take a minute.)")
    start_time = time.perf_counter()

    for i in range(iter_):
        build_atomistic_directed(polymers[i % len(polymers)])

    run_time = time.perf_counter() - start_time
    print(f"Done graph time test:{datetime.datetime.now().time()}")
    return run_time/iter_ * 1000  # micro-seconds


def performance_test():
    polymer_string = [
        "CC{[>][<]CC(C)[>][<]}CC(C)=C",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        "CC{[>][$]CC[$],[$]CC(CC)[$][<]}O",
        "CC{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}O",
        "CC{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}F",
        "C{[$][$]CC[$],[$]CC(CC)[$][$]}CC",
        "O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}CC"
    ]
    time_iter = 2_000
    memory_iter = 1000
    time_out= time_bigsmiles_parsing_graph(polymer_string, time_iter)
    print(time_out)


def main():
    polymer = "O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCC(C)N[>][<]}CF"  # wrong end group can terminate on other
    # polymer = "OC{[>][<]CC(C{[>][<]CCO[>][<]}CN)[>][<]}CC"
    # polymer = "CCC(C){[$][$]CC(C1CCCCC1)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][$]}[H]"
    # polymer = "C{[$][$]CC[$],[$]CC(CC[$])[$][$]}O"
    # polymer ="CC(CC){[<1][>1]CC(C)[<2][>2]}CCO"  # Fails
    # polymer = "C{[<][>]CC[>2],[<2]N[>][<]}O"   # Fails
    # polymer = "C{[<]CC(C[>])C(C[<])CO[>]}O"
    # polymer = "F{[<][>]CO[>],[>]C(N[<])C[<][>]}F"  # end group wrong can terminate on other
    polymer = "F{[<][>]OO[>],[>]C(N[<])C[<][>]}S"  # end group wrong
    # polymer = "{[][>]C([>])([>]),[<]OO[>][>]}CB"  # error
    # polymer = "F{[<1][>1]C([<2])C[<1],[>2]OO[<2][>2]}N"  #
    # polymer = "F{[<][>]C([<])C[<],[>]OO[<][>]}N"  #
    # polymer = "O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCC(C)N[>][<]}CF"
    # polymer ="F{[<][>]CC(F)[<],[>]CCO[<][>]}P"

    # polymer = "CC{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}O"


    graph = build_atomistic_directed(polymer)
    print(graph)
    draw(graph)


if __name__ == "__main__":
    main()
    # performance_test()