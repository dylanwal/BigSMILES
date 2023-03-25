import math

from bigsmiles.distributions.base import Distribution


def _plotting_general(x, y, name: str, log_scale: bool = True, normalize: bool = True, fig=None):
    try:
        import plotly.graph_objs as go
    except ImportError:
        raise ImportError("'plotly' is an optional package and needs to be installed. "
                          "\n'pip install plotly'")

    if fig is None:  # create figure is one is not provided
        fig = go.Figure()

    if normalize:
        y = y / max(y)

    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name=str(name)))

    if log_scale:
        fig.update_xaxes(type="log")

    return fig


def _add_layout(fig):
    fig.layout.legend.font.size = 10
    fig.update_layout(autosize=False, width=800, height=500, font=dict(family="Arial", size=18, color="black"),
                      plot_bgcolor="white", showlegend=True, legend=dict(x=.02, y=.95))
    fig.update_xaxes(tickprefix="<b>", ticksuffix="</b>", showline=True, linewidth=5, mirror=True, linecolor='black',
                     ticks="outside", tickwidth=4, showgrid=False, gridwidth=1, gridcolor="lightgray")
    fig.update_yaxes(tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray")


def plot_x_i(dis: Distribution, log_scale: bool = True, normalize: bool = True, fig=None):
    """
    for quick visualization of distribution x_i

    Parameters
    ----------
    dis:
        Distribution
    log_scale:
        True: x-axis log-scale; False: x-axis is linear scale
    normalize:
        True: max set 1
    fig:
        plotly figure you want to add the trace to

    Returns
    -------
    fig:
        plotly figure
    """
    fig_out = _plotting_general(dis.mw_i, dis.x_i(dis.mw_i), str(dis), log_scale, normalize, fig)

    if fig is None:
        _add_layout(fig_out)
        fig_out.layout.xaxis.title = "<b>molecular weight, <i>mw<sub>i</sub></i> (g/mol)</b>"
        if normalize:
            fig_out.layout.yaxis.title = "<b>normalized mole fraction, <i>x<sub>i</sub></i></b>"
        else:
            fig_out.layout.yaxis.title = "<b>mole fraction, <i>x<sub>i</sub></i></b>"

    if hasattr(dis, "_plot_range"):
        fig_out.update_xaxes(range=[math.log10(i) for i in dis._plot_range])

    return fig_out


def plot_x_i_pmd(dis: Distribution, log_scale: bool = True, normalize: bool = True, fig=None):
    """
    for quick visualization of distribution x_i_pmd

    Parameters
    ----------
    dis:
        Distribution
    log_scale:
        True: x-axis log-scale; False: x-axis is linear scale
    normalize:
        True: max set 1
    fig:
        plotly figure you want to add the trace to

    Returns
    -------
    fig:
        plotly figure
    """
    fig_out = _plotting_general(dis.N_i, dis.x_i_pmd(), str(dis), log_scale, normalize, fig)

    if fig is None:
        _add_layout(fig_out)
        fig_out.layout.xaxis.title = "<b>chain length, <i>N<sub>i</sub></i></b>"
        if normalize:
            fig_out.layout.yaxis.title = "<b>normalized mole fraction, <i>x<sub>i</sub></i></b>"
        else:
            fig_out.layout.yaxis.title = "<b>mole fraction, <i>x<sub>i</sub></i></b>"

    if hasattr(dis, "_plot_range"):
        fig_out.update_xaxes(range=[math.log10(i) for i in dis._plot_range])

    return fig_out


def plot_w_i(dis: Distribution, log_scale: bool = True, normalize: bool = True, fig=None):
    """
    for quick visualization of distribution w_i

    Parameters
    ----------
    dis:
        Distribution
    log_scale:
        True: x-axis log-scale; False: x-axis is linear scale
    normalize:
        True: max set 1
    fig:
        plotly figure you want to add the trace to

    Returns
    -------
    fig:
        plotly figure
    """
    fig_out = _plotting_general(dis.mw_i, dis.w_i(dis.mw_i), str(dis), log_scale, normalize, fig)

    if fig is None:
        _add_layout(fig_out)
        fig_out.layout.xaxis.title = "<b>molecular weight, <i>mw<sub>i</sub></i> (g/mol)</b>"
        if normalize:
            fig_out.layout.yaxis.title = "<b>normalized weight fraction, <i>w<sub>i</sub></i></b>"
        else:
            fig_out.layout.yaxis.title = "<b>weight fraction, <i>w<sub>i</sub></i></b>"

    if hasattr(dis, "_plot_range"):
        fig_out.update_xaxes(range=[math.log10(i) for i in dis._plot_range])

    return fig_out


def plot_x_i_cdf(dis: Distribution, fig=None):
    """
    for quick visualization of distribution x_i_cdf

    Parameters
    ----------
    dis:
        Distribution
    fig:
        plotly figure you want to add the trace to

    Returns
    -------
    fig:
        plotly figure
    """
    x, y = dis.x_i_cdf()
    fig_out = _plotting_general(x, y, str(dis), False, True, fig)

    if fig is None:
        _add_layout(fig_out)
        fig_out.layout.xaxis.title = "<b>molecular weight, <i>mw<sub>i</sub></i> (g/mol)</b>"
        fig_out.layout.yaxis.title = "<b> cumulative mole fraction </b>"

    if hasattr(dis, "_plot_range"):
        fig_out.update_xaxes(range=[0, 1])

    return fig_out


def plot_w_i_cdf(dis: Distribution, fig=None):
    """
    for quick visualization of distribution w_i_cdf

    Parameters
    ----------
    dis:
        Distribution
    fig:
        plotly figure you want to add the trace to

    Returns
    -------
    fig:
        plotly figure
    """
    x, y = dis.w_i_cdf()
    fig_out = _plotting_general(x, y, str(dis), False, True, fig)

    if fig is None:
        _add_layout(fig_out)
        fig_out.layout.xaxis.title = "<b>molecular weight, <i>mw<sub>i</sub></i> (g/mol)</b>"
        fig_out.layout.yaxis.title = "<b> cumulative wight fraction </b>"

    if hasattr(dis, "_plot_range"):
        fig_out.update_xaxes(range=[0, 1])

    return fig_out
