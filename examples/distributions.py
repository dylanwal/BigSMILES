import bigsmiles

dis = bigsmiles.distributions.LogNormal(Mn=10_000, D=1.12, repeat_MW=104.15)
# dis = bigsmiles.distributions.Poisson(N=100, repeat_MW=104.12)

print('Mn: ', dis.Mn)
print('Mw: ', dis.Mw)
print('D: ', dis.D)
print('skew in mw: ', dis.skew_mw)
print('standard deviation in mw: ', dis.std_mw)
print('kurtosis in mw: ', dis.kurtosis_mw)
print('peak mw: ', dis.peak_mw)
print('chain length', dis.N)

fig = bigsmiles.distributions.plot_w_i(dis)
fig.show()
# fig.write_image("quickstart_distributions.svg")


# Picking mw
import numpy as np
import plotly.graph_objs as go

n = 1000
mw = dis.draw_mw(n)

# histogram
bins = np.logspace(np.log10(np.min(mw)), np.log10(np.max(mw)), 20)
y, edges = np.histogram(mw, bins=bins)
# normalize y
y = y / np.max(y)
# center point of x
x = (edges[:-1] + edges[1:]) / 2
width = edges[1:] - edges[:-1]

fig = bigsmiles.distributions.plot_w_i(dis)
fig.add_trace(
    go.Bar(
        x=x,
        y=y,
        width=width,
        name="drawn mw"
    )
)

fig.show()
