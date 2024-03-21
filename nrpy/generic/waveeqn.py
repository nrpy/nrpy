"""
The waveequation! It can't be solved too many times.
"""

from typing import cast
from sympy import IndexedBase, Idx, symbols, Basic, Expr, Indexed, Symbol, Function
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
gf.fill_in(p_t[-i], lambda _,i: mkSymbol(f"pD{i}"))
gf.fill_in(u_d[-i])
gf.fill_in(p_d[-i,-j], lambda _,i,j: mkSymbol(f"pD{i}_dD{j}"))


eqnlist = EqnList()
for ii in range(3):
    term = do_subs(u_d[-ii-1], gf.subs)
    eqnlist.add_input(cast(Symbol, term))

eqnlist.add_input(mkSymbol("p"))
div2 = Function("div2")
def mkdiv2(var, a, b):
    return div2(mkSymbol(var), sympify(a), sympify(b))
eqnlist.add_eqn(mkSymbol("pD0_dD0"), mkdiv2("p", 0, 0))
eqnlist.add_eqn(mkSymbol("pD1_dD1"), mkdiv2("p", 1, 1))
eqnlist.add_eqn(mkSymbol("pD2_dD2"), mkdiv2("p", 2, 2))

eqn1 = mkEq(p_t[-j], u_d[-j])
for eqn in gf.expand_eqn(eqn1):
    eqnlist.add_output(eqn.lhs)
    eqnlist.add_eqn(eqn.lhs, eqn.rhs)
eqnlist.add_output(cast(Symbol, u_t))
eqn2 = mkEq(u_t, gf.expand(g[i,j]*p_d[-i,-j]))

eqnlist.add_eqn(eqn2.lhs, eqn2.rhs)

eqnlist.diagnose()
eqnlist.dump()
eqnlist.cse()
eqnlist.diagnose()
eqnlist.dump()

def generate_cactus_thorn(
    eqnlist : EqnList
    project_dir : str,
    thorn_name : str,
    inherits : Optional[str]=None)->None:
    construct_interface_ccl(
        project_dir=project_dir,
        thorn_name=thorn_name,
        inherits=inherits,
        USES_INCLUDEs="", # Not sure
        is_evol_thorn=True,
        enable_NewRad=False)

generate_cactus_thorn(eqnlist, "wavetoy", "WaveToyProj")
