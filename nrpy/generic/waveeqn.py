"""
The waveequation! It can't be solved too many times.
"""

from sympy import IndexedBase, Idx, symbols
from nrpy.generic.sympywrap import *
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.use_indices import expand, add_sym, fill_in, GF

i, j, k = mkIdxs('i j k')

g = mkIndexedBase('g',(3,3))
add_sym(g[i,j],i,j)

p = mkIndexedBase('p',(3,))
u_d = mkIndexedBase('u_d',(3,))

fill_in(g[i,j], lambda i,j: 1 if i==j else 0)
fill_in(p[-i], lambda i: "p_"+["x","y","z"][i])
fill_in(u_d[-i], lambda i: "u_dD"+str(i))

dvdt = g[i,j]*p[-i]*u_d[-j]
print(dvdt)
res = expand(dvdt)
print(res)

#dudt = dpdt

gf = GF()
h = gf.decl("h",[-i,-j])
print(h)
