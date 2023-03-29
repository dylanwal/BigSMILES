
try:
    import networkx as nx
except ImportError:
    raise ImportError("To use this feature install networkx. (pip install networkx)")
try:
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("To use this feature install matplotlib. (pip install matplotlib)")


colors = {
    "C": (0,0,0),
    "O": (199/255, 41/255,  24/255),
    "N": (8/255, 0/255, 128/255),
    "F": (8/255, 138/255, 28/255),
    "S": (205/255, 212/255, 8/255),
    "Q": (137/255, 8/255, 212/255),
    "<": (137/255, 8/255, 212/255),
    ">": (137/255, 8/255, 212/255)
}


def get_atom_colors(graph) -> list:
    atom_colors = []

    for node in graph.nodes:
        if node[0] in colors:
            atom_colors.append(colors[node[0]])
        else:
            atom_colors.append(colors["C"])

    return atom_colors


edge_color_map = {
    "": (0,0,0),
    "=": (1,0,0),
    "#": (0,1,0)
}


def draw(graph):
    # edge_color = [edge_color_map[edge['symbol']] for _, _, edge in graph.edges(data=True)]

    plt.close()
    plt.ylim([-1.1, 1.1])
    plt.xlim([-1.1, 1.1])
    plt.axis('off')
    nx.draw_networkx(graph,
                     pos=nx.kamada_kawai_layout(graph, pos=nx.random_layout(graph)),
                     # pos=kamada_kawai_layout(graph),
                     # pos=nx.spring_layout(graph, iterations=200),
                     # pos = nx.circular_layout(graph),
                     # edge_color=edge_color,
                     node_size=250,
                     arrowsize=10,
                     arrows=True,
                     node_color=get_atom_colors(graph),
                     width=2,
                     with_labels=True,
                     font_size=6,
                     font_color="w",
                     font_family="Arial Rounded MT Bold"
                     )
    plt.show()


def draw_plotly(graph):
    import plotly.graph_objects as go

    pos = nx.kamada_kawai_layout(graph, pos=nx.random_layout(graph))
    edge_x = []
    edge_y = []
    for edge in graph.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='text',
        hovertext=['<br>'.join([f"{k}: {v}" for k, v in data.items()]) for _, _, data in graph.edges(data=True)],
        mode='lines')

    node_x = []
    node_y = []
    for node in graph.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        # hovertemplate='%{text}',
        hovertext= ['<br>'.join([f"{k}: {v}" for k, v in data.items()]) for _, data in graph.nodes(data=True)],
        marker=dict(
            color=[],
            size=10,
            line_width=2))

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        showlegend=False,
                        hovermode='closest',
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.write_html('temp.html', auto_open=True)
