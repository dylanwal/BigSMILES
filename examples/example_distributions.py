import plotly.graph_objs as go

import bigsmiles.reference_data.distributions as distributions


dis = distributions.LogNormal(10_000, 1.2)
fig = dis.plot()
fig.write_html('temp.html', auto_open=True)




