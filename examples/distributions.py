import bigsmiles.distributions as distributions


def main():
    repeat_MW = 104.15  # styrene
    Mn = 10_000
    D = 1.1
    dis = distributions.Uniform(Mn, D, repeat_MW)
    fig = distributions.plot_x_i(dis)
    fig = distributions.plot_w_i(dis, fig=fig)
    fig = distributions.plot_x_i_pmd(dis, fig=fig)
    fig.show()

    print(dis.details())


if __name__ == "__main__":
    main()
