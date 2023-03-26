
import bigsmiles
import bigsmiles.distributions as big_dist

polymer = bigsmiles.Polymer("[H]{[>][<]CC(C1=CC=CC=C1)[>][<]}[H]")
dis = big_dist.LogNormal(Mn=10_000, D=1.12)

polymer.add_distribution(polymer.stochastic_objects, dis)
