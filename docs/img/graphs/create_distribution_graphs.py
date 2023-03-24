"""

create distribution graphs

"""
import math

import plotly.graph_objs as go
import numpy as np

import bigsmiles.distributions as distributions


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
    # fig.layout.legend.font.size = 10
    # fig.update_layout(autosize=False, width=800, height=500, font=dict(family="Arial", size=18, color="black"),
    #                   plot_bgcolor="white", showlegend=True, legend=dict(x=.02, y=.95))
    # fig.update_xaxes( tickprefix="<b>", ticksuffix="</b>", showline=True,
    #                  linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
    #                  gridwidth=1, gridcolor="lightgray")
    # fig.update_yaxes(tickprefix="<b>", ticksuffix="</b>", showline=True,
    #                  linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
    #                  gridwidth=1, gridcolor="lightgray")


def main():
    n = 6
    repeat_MW = 104.15  # styrene
    Mn = np.logspace(3, 6, n)
    D = np.linspace(1.03, 1.8, n)

    # LogNormal
    for i in range(n):
        dis = distributions.LogNormal(Mn[i], D[i])
        if i == 0:
            fig = distributions.plot_w_i(dis)
        else:
            fig = distributions.plot_w_i(dis, fig=fig)

    format_figure(fig)
    fig.write_image("LogNormal.svg")

    # SchulzZimm
    for i in range(n):
        dis = distributions.SchulzZimm(Mn[i], D[i])
        if i == 0:
            fig = distributions.plot_w_i(dis)
        else:
            fig = distributions.plot_w_i(dis, fig=fig)

    format_figure(fig)
    fig.write_image("SchulzZimm.svg")

    # gauss
    d_ = np.linspace(1.01, 1.2, n)
    for i in range(n):
        dis = distributions.Gaussian(Mn[i], d_[i])
        if i == 0:
            fig = distributions.plot_w_i(dis)
        else:
            fig = distributions.plot_w_i(dis, fig=fig)

    format_figure(fig)
    fig.write_image("Gaussian.svg")

    # uniform
    wide = np.logspace(3, 6, n)
    for i in range(n):
        dis = distributions.Uniform(500, wide[i])
        if i == 0:
            fig = distributions.plot_x_i(dis)
        else:
            fig = distributions.plot_x_i(dis, fig=fig)

    format_figure(fig)
    fig.write_image("Uniform.svg")

    # CustomDistribution
    data = np.genfromtxt("../../../examples/data/GPC_data.csv", delimiter=',')
    for i in range(1, len(data[0])):
        x = data[:, 0]
        y = data[:, i]
        dis = distributions.CustomDistribution(mw_i=x, w_i=y)
        if i == 1:
            fig = distributions.plot_w_i(dis)
        else:
            fig = distributions.plot_w_i(dis, fig=fig)

    format_figure(fig)
    fig.write_image("CustomDistribution.svg")

    # FlorySchulz
    conversion = [0.40, 0.95, 0.99, 0.995, 0.999, 0.9995]
    for i in range(len(conversion)):
        dis = distributions.FlorySchulz(conversion[i], repeat_MW)
        if i == 0:
            fig = distributions.plot_x_i_pmd(dis)
        else:
            fig = distributions.plot_x_i_pmd(dis, fig=fig)

    format_figure(fig)
    fig.write_image("FlorySchulz.svg")

    # Poisson
    N = np.linspace(20, 800, n, dtype='int')
    for i in range(n):
        dis = distributions.Poisson(N[i], repeat_MW)
        if i == 0:
            fig = distributions.plot_x_i_pmd(dis)
        else:
            fig = distributions.plot_x_i_pmd(dis, fig=fig)

    format_figure(fig)
    fig.write_image("Poisson.svg")

    print("done")


if __name__ == "__main__":
    main()
