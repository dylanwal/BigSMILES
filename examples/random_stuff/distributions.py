import numpy as np
import plotly.graph_objs as go


import bigsmiles.distributions as distributions


def main():
    repeat_MW = 104.15  # styrene
    Mn = 10_000
    D = 1.1
    # dis = distributions.Uniform(Mn, D, repeat_MW)
    # dis = distributions.Uniform(2000, 20000, repeat_MW)
    fig = distributions.plot_x_i(dis)
    fig = distributions.plot_w_i(dis, fig=fig)
    fig = distributions.plot_x_i_pmd(dis, fig=fig)
    fig.write_html('temp.html', auto_open=True)

    print(dis.details())


def main2():
    n = 6
    repeat_MW = 104.15  # styrene
    Mn = np.logspace(3, 6, n)
    D = np.linspace(1.03, 1.8, n)

    conversion = [0.40, 0.95, 0.99, 0.995, 0.999, 0.9995]
    for i in range(len(conversion)):
        dis = distributions.FlorySchulz(conversion[i], repeat_MW)
        dis.N_i
        if i == 0:
            fig = distributions.plot_x_i_pmd(dis)
        else:
            fig = distributions.plot_x_i_pmd(dis, fig=fig)

    fig.write_html("temp.html", auto_open=True)


if __name__ == "__main__":
    main2()
