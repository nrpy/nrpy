"""
The waveequation! It can't be solved too many times.
"""

from typing import cast, Optional, Set, Any
from sympy import IndexedBase, Idx, symbols, Basic, Expr, Indexed, Symbol, Function
from nrpy.generic.sympywrap import *
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.use_indices import *
from nrpy.infrastructures.ETLegacy.interface_ccl import construct_interface_ccl
from nrpy.helpers.colorize_text import colorize

def flat_metric(out:Expr, ni:Idx, nj:Idx)->Expr:
    i = to_num(ni)
    j = to_num(nj)
    if i==j:
        return sympify(1)
    else:
        return sympify(0)

# Create a set of grid functions
gf = GF()

# Declare gfs
p   = gf.decl("p", [li], "VVC")
p_t = gf.decl("p_t", [li], "VVC")
p_d = gf.decl("p_d",[li,lj], "VVC")
u   = gf.decl("u", [], "VVC")
u_t = gf.declscalar("u_t")
u_d = gf.decl("u_d", [ui], "VVC")

siter2 = gf.decl("siter2", [li,lj])
gf.add_sym(siter2[li,lj], li, lj)
iter1 = gf.decl("iter1", [li])

# Declare the metric
g = gf.decl("g",  [li,lj])
gf.add_sym(g[li,lj],li,lj)

# Declare params
spd = gf.add_param("spd", default=1.0, desc="The wave speed")

# Fill in values
gf.fill_in(g[li,lj], flat_metric)
gf.fill_in(g[ui,uj], flat_metric)

# Fill in with defaults
gf.fill_in(p[li], lambda _,i: mkSymbol(f"pD{to_num(i)}"))
gf.fill_in(p_t[li], lambda _,i: mkSymbol(f"p_tD{to_num(i)}"))


# Fill in the deriv variables with a function call
#
div1 = gf.declfun("div1", True)
divx = gf.declfun("divx", True)
divy = gf.declfun("divy", True)
divz = gf.declfun("divz", True)
gf.fill_in(p_d[li,lj], lambda _,i,j: div1(p[j],i))
gf.fill_in(u_d[li], lambda _,i: div1(u,i))

def to_div(out:Expr, j:Idx)->Expr:
    n = to_num(j)
    ret : Any
    if n==0:
        ret = divx(*out.args[:-1])
    elif n==1:
        ret = divy(*out.args[:-1])
    elif n==2:
        ret = divz(*out.args[:-1])
    else:
        assert False
    return cast(Expr, ret)

def to_div2(out:Expr, i:Idx, j:Idx)->Expr:
    n = to_num(j)
    ret : Any
    if n==0:
        ret = divx(p[i])
    elif n==1:
        ret = divy(p[i])
    elif n==2:
        ret = divz(p[i])
    else:
        assert False
    return cast(Expr, ret)

gf.fill_in(iter1[lj], alt=div1(u,lj), f=to_div)
gf.fill_in(siter2[li,lj], alt=div1(p[li],lj), f=to_div2)

res = gf.do_subs(spd*g[ui,uj]*div1(p[lj],li))

# Add the equations we want to evolve.
gf.add_eqn(p_t[lj], spd*div1(u,lj), "EVO")
gf.add_eqn(u_t, spd*g[ui,uj]*div1(p[lj],li), "EVO")

# Ensure the equations make sense
gf.diagnose()

# Display the equations in final form
gf.dump()

# Perform cse
gf.cse()

# Display again in case there are changes
gf.dump()

gf.show_tensortypes()
