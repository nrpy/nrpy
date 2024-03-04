"""
The waveequation! It can't be solved too many times.
"""

from sympy import IndexedBase, Idx, symbols
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.use_indices import expand, add_sym, fill_in, GF

i, j, k = symbols('i j k', cls=Idx)

g = IndexedBase('g')
add_sym(g[i,j],i,j)

p = IndexedBase('p')
u_d = IndexedBase('u_d')

fill_in(g[i,j], lambda i,j: 1 if i==j else 0)
fill_in(p[-i], lambda i: "p_"+["x","y","z"][i])
fill_in(u_d[-i], lambda i: "u_dD"+str(i))

dvdt = g[i,j]*p[-i]*u_d[-j]
print(dvdt)
res = expand(dvdt)
print(res)

#dudt = dpdt

gf = GF()
gf.decl("h",[-i,-j])
print(h)
