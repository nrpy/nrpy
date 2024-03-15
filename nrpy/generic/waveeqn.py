"""
The waveequation! It can't be solved too many times.
"""

from sympy import IndexedBase, Idx, symbols, Basic, Expr, Indexed
from nrpy.generic.sympywrap import *
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.use_indices import expand, add_sym, fill_in, GF, expand_free_indices

i, j, k = mkIdxs('i j k')

def flat_metric(out:Indexed, i:int, j:int)->Expr:
    if i==j:
        return sympify(1)
    else:
        return sympify(0)

# Create a set of grid functions
gf = GF()

# Declare stuff
p_t = gf.decl("p_t",[i])
p_d = gf.decl("p_d",[i,j])
u_t = gf.decl("u_t",[])
u_d = gf.decl("u_d",[i])
g   = gf.decl("g",  [i,j])

# Fill in values
gf.fill_in(g[i,j], flat_metric)
# Fill in with defaults
gf.fill_in(p_t[-i])
gf.fill_in(u_d[-i])
gf.fill_in(p_d[-i,-j])


eqnlist = EqnList()
for ii in range(3):
    term = do_subs(u_d[-ii-1], gf.subs)
    eqnlist.add_input(term)
    print(term)

eqn1 = mkEq(p_t[-j], u_d[-j])
for eqn in gf.expand_eqn(eqn1):
    print(eqn)
    eqnlist.add_output(eqn.lhs)
    eqnlist.add_eqn(eqn.lhs, eqn.rhs)
eqnlist.add_output(u_t)
eqn2 = mkEq(u_t, gf.expand(g[i,j]*p_d[-i,-j]))

print(eqn2)

eqnlist.add_eqn(eqn2.lhs, eqn2.rhs)

eqnlist.diagnose()
print(eqnlist.inputs)
print(eqnlist.outputs)
