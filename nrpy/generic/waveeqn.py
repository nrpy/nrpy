"""
The waveequation! It can't be solved too many times.
"""

from typing import cast, Optional
from sympy import IndexedBase, Idx, symbols, Basic, Expr, Indexed, Symbol, Function
from nrpy.generic.sympywrap import *
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.use_indices import expand, add_sym, fill_in, GF, expand_free_indices
from nrpy.infrastructures.ETLegacy.interface_ccl import construct_interface_ccl

i, j, k = mkIdxs('i j k')

def flat_metric(out:Indexed, i:int, j:int)->Expr:
    if i==j:
        return sympify(1)
    else:
        return sympify(0)

# Create a set of grid functions
gf = GF()

# Declare gfs
p   = gf.decl("p",[i])
p_t = gf.decl("p_t",[i])
p_d = gf.decl("p_d",[i,j])
u   = gf.decl("u",[])
u_t = gf.decl("u_t",[])
u_d = gf.decl("u_d",[i])

# Declare the metric
g   = gf.decl("g",  [i,j])

# Declare params
spd = gf.add_param("spd", default=1.0, desc="The wave speed")

# Fill in values
gf.fill_in(g[i,j], flat_metric)

# Fill in with defaults
gf.fill_in(p_t[-i], lambda _,i: mkSymbol(f"p_tD{i}"))
gf.fill_in(p[-i], lambda _,i: mkSymbol(f"pD{i}"))


# Fill in the deriv variables with a function call
#
div2 = Function("div2")
def mkdiv2(var:str, a:int, b:int)->Symbol:
    return cast(Symbol, div2(mkSymbol(var), sympify(a), sympify(b)))
gf.fill_in(p_d[-i,-j], lambda _,i,j: mkdiv2("p",i,j))
#
div1 = Function("div1")
def mkdiv1(var:str, a:int)->Symbol:
    return cast(Symbol, div1(mkSymbol(var), sympify(a)))
gf.fill_in(u_d[-i], lambda _,i: mkdiv1("u",i))

# Add the equations we want to evolve.
gf.add_eqn(p_t[-j], spd*u_d[-j])
gf.add_eqn(u_t, spd*g[i,j]*p_d[-i,-j] )

# Ensure the equations make sense
gf.diagnose()

# Display the equations in final form
gf.dump()

# Perform cse
gf.cse()

# Display again in case there are changes
gf.dump()

def generate_cactus_thorn(
    gf : GF,
    project_dir : str,
    thorn_name : str,
    inherits : Optional[str]=None)->None:
    for gfn in gf.gfs:
        if gfn in gf.eqnlist.outputs or gfn in gf.eqnlist.inputs:
            print(colorize(gfn,"magenta"))
    #construct_interface_ccl(
    #    project_dir=project_dir,
    #    thorn_name=thorn_name,
    #    inherits=inherits,
    #    USES_INCLUDEs="", # Not sure
    #    is_evol_thorn=True,
    #    enable_NewRad=False)

generate_cactus_thorn(gf, "wavetoy", "WaveToyProj")
