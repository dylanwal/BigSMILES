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
    fig = _plotting_general(dis.mw_i, dis.x_i(dis.mw_i), str(dis), log_scale, normalize, fig)

    if hasattr(dis, "_plot_range"):
        fig.update_xaxes(range=[math.log10(i) for i in dis._plot_range])

    return fig


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
    fig = _plotting_general(dis.N_i, dis.x_i_pmd(), str(dis), log_scale, normalize, fig)

    if hasattr(dis, "_plot_range"):
        fig.update_xaxes(range=[math.log10(i) for i in dis._plot_range])

    return fig


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
    fig = _plotting_general(dis.mw_i, dis.w_i(dis.mw_i), str(dis), log_scale, normalize, fig)

    if hasattr(dis, "_plot_range"):
        fig.update_xaxes(range=[math.log10(i) for i in dis._plot_range])

    return fig


def plot_x_i_cdf(dis: Distribution, normalize: bool = True, fig=None):
    """
    for quick visualization of distribution x_i_cdf

    Parameters
    ----------
    dis:
        Distribution
    normalize:
        True: max set 1
    fig:
        plotly figure you want to add the trace to

    Returns
    -------
    fig:
        plotly figure
    """
    x, y = dis.x_i_cdf()
    fig = _plotting_general(x, y, str(dis), False, normalize, fig)

    if hasattr(dis, "_plot_range"):
        fig.update_xaxes(range=[0, 1])

    return fig


def plot_w_i_cdf(dis: Distribution, normalize: bool = True, fig=None):
    """
    for quick visualization of distribution w_i_cdf

    Parameters
    ----------
    dis:
        Distribution
    normalize:
        True: max set 1
    fig:
        plotly figure you want to add the trace to

    Returns
    -------
    fig:
        plotly figure
    """
    x, y = dis.w_i_cdf()
    fig = _plotting_general(x, y, str(dis), False, normalize, fig)

    if hasattr(dis, "_plot_range"):
        fig.update_xaxes(range=[0, 1])

    return fig
