"""

create distribution graphs

"""
import math

import plotly.graph_objs as go
import numpy as np

import bigsmiles.data_structures.distributions as distributions


color_list = [
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
    reset_colors(fig, color_list)
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
    D = np.linspace(1.03, 1.8, n)

    # LogNormal
    fig = go.Figure()
    labels = []
    for i in range(n):
        dis = distributions.LogNormal(Mn[i], D[i])
        fig = dis.plot_pdf(fig=fig)
        labels.append(str(dis))

    format_figure(fig)
    fig.write_image("LogNormal.svg")

    # SchulzZimm
    fig = go.Figure()
    labels = []
    for i in range(n):
        dis = distributions.SchulzZimm(Mn[i], D[i])
        fig = dis.plot_pdf(fig=fig)
        labels.append(str(dis))

    format_figure(fig)
    fig.write_image("SchulzZimm.svg")

    # gauss
    fig = go.Figure()
    labels = []
    d_ = np.linspace(1.01, 1.2, n)
    for i in range(n):
        dis = distributions.Gauss(Mn[i], d_[i])
        fig = dis.plot_pdf(fig=fig)
        labels.append(str(dis))

    format_figure(fig)
    fig.write_image("Gauss.svg")

    # uniform
    fig = go.Figure()
    labels = []
    wide = np.logspace(3, 6, n)
    for i in range(n):
        dis = distributions.Uniform(500, wide[i])
        fig = dis.plot_pdf(fig=fig)
        labels.append(str(dis))

    format_figure(fig)
    fig.write_image("Uniform.svg")

    # CustomDistribution
    fig = go.Figure()
    labels = []
    data = np.genfromtxt("../../../examples/data/GPC_data.csv", delimiter=',')
    for i in range(1, len(data[0])):
        x = data[:, 0]
        y = data[:, i]
        dis = distributions.CustomDistribution(mw_i=x, pdf=y)
        fig = dis.plot_pdf(fig=fig)
        labels.append(str(dis))

    format_figure(fig)
    fig.update_xaxes(range=(math.log10(300), math.log10(500_000)))
    fig.write_image("CustomDistribution.svg")

    # FlorySchulz
    fig = go.Figure()
    labels = []
    conversion = [0.40, 0.95, 0.99, 0.995, 0.999, 0.9995]
    for conv in conversion:
        dis = distributions.FlorySchulz(conv, repeat_MW)
        fig = dis.plot_pmf(fig=fig)
        labels.append(str(dis))

    format_figure(fig)
    fig.update_xaxes(range=(0, 4))
    fig.layout.xaxis.title = "<b>chain length</b>"
    fig.write_image("FlorySchulz.svg")

    # Poisson
    fig = go.Figure()
    labels = []
    N = np.linspace(20, 800, n, dtype='int')
    for i in range(n):
        dis = distributions.Poisson(N[i], repeat_MW)
        fig = dis.plot_pmf(fig=fig)
        labels.append(str(dis))

    format_figure(fig)
    fig.update_xaxes(range=(0, 4))
    fig.layout.xaxis.title = "<b>chain length</b>"
    fig.write_image("Poisson.svg")

    print("done")


if __name__ == "__main__":
    main()
