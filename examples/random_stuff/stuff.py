
import bigsmiles.distributions as distributions

dis = distributions.FlorySchulz(conversion=0.999, repeat_MW=104)

print(dis.N)
print(dis.std_mw/dis.repeat_MW)
print(dis.skew_mw/dis.repeat_MW)
