"""
The waveequation! It can't be solved too many times.
"""

from typing import cast, Optional, Set
from sympy import IndexedBase, Idx, symbols, Basic, Expr, Indexed, Symbol, Function
from nrpy.generic.sympywrap import *
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.use_indices import expand, add_sym, fill_in, GF, expand_free_indices
from nrpy.infrastructures.ETLegacy.interface_ccl import construct_interface_ccl
from nrpy.helpers.colorize_text import colorize

i, j, k = mkIdxs('i j k')

def flat_metric(out:Indexed, i:int, j:int)->Expr:
    if i==j:
        return sympify(1)
    else:
        return sympify(0)

# Create a set of grid functions
gf = GF()
m = True
gf.do_div1repl = m

# Declare gfs
p   = gf.decl("p",[-i])
p_t = gf.decl("p_t",[-i])
p_d = gf.decl("p_d",[-i,-j])
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
gf.fill_in(p[i], lambda _,i: mkSymbol(f"pD{i}"))
gf.fill_in(p_t[-i], lambda _,i: mkSymbol(f"p_tD{i}"))


# Fill in the deriv variables with a function call
#
div1 = Function("div1")
def mkdiv1(var:IndexedBase, a:int)->Symbol:
    return cast(Symbol, div1(var, sympify(a)))

if m:
    gf.fill_in(p_d[-i,-j], lambda _,i,j: mkdiv1(p[j],i))
    gf.fill_in(u_d[-i], lambda _,i: mkdiv1(u,i))
else:
    # div1(p[j], i) -> p_D0_dD1
    gf.fill_in(p_d[-i,-j], lambda _,i,j: div1(p[j],i), base_zero=not m)
    gf.fill_in(u_d[-i], lambda _,i: div1(u,i), base_zero=not m)

# Add the equations we want to evolve.
gf.add_eqn(p_t[-j], spd*div1(u,-j), "EVO")
gf.add_eqn(u_t, spd*g[i,j]*div1(p[-j],-i), "EVO")

# Ensure the equations make sense
gf.diagnose()

# Display the equations in final form
gf.dump()

# Perform cse
gf.cse()

# Display again in case there are changes
gf.dump()

# Probably want to have multiple GF objects
#gf.set_function("wave_evol")
#gf.set_schedule("MoL_CalcRHS")

def generate_cactus_thorn(
    gf : GF,
    project_dir : str,
    thorn_name : str,
    inherits : Optional[str]=None)->None:

    output : Set[str] = set()
    for gfn in gf.gfs:
        gfs = mkSymbol(gfn)
        if gfs in gf.eqnlist.outputs or gfs in gf.eqnlist.inputs:
            if gfn in gf.base_of:
                gof = gf.base_of[gfn]
                if gof in output:
                    continue
                output.add(gof)
                #print(f"base of '{gfn}' is '{gof}'")
                #print(f">> {gof} {{ {gf.groups[gof]} }}")
                gplist = list(gf.groups[gof])
                gflist = " ".join(gplist)
                comment = gf.defn[gof]
                centering = ""
                checkpointing = ""
                tags = ""
                props = gf.props[gfn]
                if len(props) == 2:
                    if props[0] < 0 and props[1] < 0:
                        if len(gplist) == 6: # and dim == 3
                            tags += ' tensortypealias="DD_sym"'
                        else:
                            tags += ' tensortypealias="DD"'
                    elif props[0] > 0 and props[1] > 0:
                        if len(gplist) == 6: # and dim == 3
                            tags += ' tensortypealias="UU_sym"'
                        else:
                            tags += ' tensortypealias="UU"'
                    elif props[0] > 0 and props[1] < 0:
                        tags += ' tensortypealias="DU"'
                    elif props[0] < 0 and props[1] > 0:
                        tags += ' tensortypealias="UD"'
                    else:
                        assert False
                elif len(props) == 1:
                    if props[0] > 0:
                        tags += ' tensortypealias="U"'
                    elif props[0] < 0:
                        tags += ' tensortypealias="D"'
                    else:
                        assert False
                elif len(props) == 0:
                    tags += ' tensortypealias="Scalar"'
                if len(tags) > 0:
                    tags=f"TAGS='{tags.strip()}'"
                print(f"CCTK_REAL {gof} TYPE=GF {centering} {tags} {{ {gflist} }} \"{comment}\"")
            else:
                tags = ""
                centering=""
                comment = gf.defn[gfn]
                tags += ' tensortypealias="Scalar"'
                if len(tags) > 0:
                    tags=f"TAGS='{tags.strip()}'"
                print(f"CCTK_REAL {gfn} TYPE=GF {centering} {tags} {{ {gfn} }} \"{comment}\"")
    #construct_interface_ccl(
    #    project_dir=project_dir,
    #    thorn_name=thorn_name,
    #    inherits=inherits,
    #    USES_INCLUDEs="", # Not sure
    #    is_evol_thorn=True,
    #    enable_NewRad=False)

generate_Carpet_thorn(gf, "wavetoy", "WaveToyProj")
generate_CarpetX_thorn(gf, "wavetoy", "WaveToyProj")
generate_StandaloneRust(gf)
