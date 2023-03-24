

import bigsmiles.distributions


dis = bigsmiles.distributions.LogNormal(Mn=10_000, D=1.12, repeat_MW=104.15)

print('Mn: ', dis.Mn)
print('Mw: ', dis.Mw)
print('D: ', dis.D)
print('skew in mw: ', dis.skew_mw)
print('standard deviation in mw: ', dis.std_mw)
print('kurtosis in mw: ', dis.kurtosis_mw)
print('peak mw: ', dis.peak_mw)
print('chain length', dis.N)

fig = bigsmiles.distributions.plot_w_i(dis)
# fig.show()
fig.write_image("quickstart_distributions.svg")
