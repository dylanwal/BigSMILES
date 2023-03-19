import plotly.graph_objs as go
import numpy as np

import bigsmiles.reference_data.distributions as distributions


color_list = [
    '#a31f34',
    '#b53133',
    '#c64230',
    '#d5552c',
    '#e36827',
    '#ef7c1f',
    '#f89113',
    '#ffa600'
    ]

color_list2 = [
    "#a31f34",
    "#ac464c",
    "#b46566",
    "#b88180",
    "#ba9d9c",
    "#b9b9b9"
]


def reset_colors(fig: go.Figure, colors: list[str, ...] | tuple[str, ...]):
    for i, plot in enumerate(fig.data):
        if isinstance(plot, go.Scatter):
            plot.line.color = colors[i]


def reset_legend_labels(fig: go.Figure, labels: list[str, ...] | tuple[str, ...]):
    for i, plot in enumerate(fig.data):
        if isinstance(plot, go.Scatter):
            plot.line.color = labels[i]


def format_figure(fig):
    reset_colors(fig, color_list2)
    # reset_legend_labels(fig, labels)
    fig.layout.legend.font.size = 10
    fig.update_layout(autosize=False, width=800, height=500, font=dict(family="Arial", size=18, color="black"),
                      plot_bgcolor="white", showlegend=True, legend=dict(x=.02, y=.95))
    fig.update_xaxes(title="<b>molecular weight (g/mol)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray", range=[1, 7])
    fig.update_yaxes(title="<b>normalized mole fraction</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray", range=[-0.05, 1.1])


def main():
    n = 6
    repeat_MW = 104.15  # styrene
    Mn = np.logspace(3, 6, n)
    D = np.linspace(1.03, 2, n)
    N = Mn/repeat_MW


    # # LogNormal
    # fig = go.Figure()
    # labels = []
    # for i in range(n):
    #     dis = distributions.LogNormal(Mn[i], D[i])
    #     fig = dis.plot_pdf(fig=fig)
    #     labels.append(str(dis))
    #
    # format_figure(fig)
    # fig.write_image("LogNormal.svg")
    #
    #
    # # FlorySchulz
    # fig = go.Figure()
    # labels = []
    # conversion = [0.8, 0.93, 0.99, 0.995, 0.999, 0.9999]
    # for conv in conversion:
    #     dis = distributions.FlorySchulz(conv, repeat_MW)
    #     fig = dis.plot_pdf(fig=fig)
    #     labels.append(str(dis))
    #
    # format_figure(fig)
    # fig.write_image("FlorySchulz.svg")


    # # SchulzZimm
    # fig = go.Figure()
    # labels = []
    # for i in range(n):
    #     dis = distributions.SchulzZimm(Mn[i], 1.1)
    #     print(dis.D)
    #     fig = dis.plot_pdf(fig=fig)
    #     labels.append(str(dis))
    #
    # format_figure(fig)
    # fig.write_image("SchulzZimm.svg")

    #
    # Poisson
    fig = go.Figure()
    labels = []
    for i in range(n):
        dis = distributions.Poisson(N[i], repeat_MW)
        fig = dis.plot_pdf(fig=fig)
        labels.append(str(dis))

    format_figure(fig)
    fig.write_image("Poisson.svg")


if __name__ == "__main__":
    main()
