import plotly.graph_objs as go

import bigsmiles.reference_data.distributions as distributions


dis = distributions.LogNormal(10_000, 1.2)
fig = dis.plot_pdf()
fig.write_html('temp.html', auto_open=True)



import numpy

def cal_Mn_D_from_xi(mw_i: np.ndarray, x_i: np.ndarray) -> tuple[float, float]:
    """ calculate Mn and D from xi vs MW data (MW goes low to high) """
    mw_n = np.sum(x_i * mw_i)
    mw_w = np.sum(x_i * mw_i**2) / mw_n
    mw_d = mw_w / mw_n
    return mw_n, mw_d



#
# from scipy import stats
# import numpy as np
# import math
#
# Mn = 10_000
# D = 1.2
# x = np.logspace(2, 6, 1000)
#
#
# def func(x_):
#     return
#
#
# std = math.sqrt(math.log(D))
# mean = math.log(Mn) - std**2/2
# distribution = stats.lognorm(s=std, loc=mean)
#
#
# fig = go.Figure()
# fig.add_trace(go.Scatter(x=x, y=distribution.pdf(x), mode='lines'))
# fig.add_trace(go.Scatter(x=x, y=func(x), mode='lines'))
# fig.update_xaxes(type="log")
# fig.show()
#
# print(np.trapz(y=func(x), x=x))
# print(np.trapz(y=distribution.pdf(x), x=x))
