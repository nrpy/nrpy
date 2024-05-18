"""
The waveequation! It can't be solved too many times.
"""

from typing import cast, Optional, Set
from sympy import IndexedBase, Idx, symbols, Basic, Expr, Indexed, Symbol, Function
from nrpy.generic.sympywrap import *
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.use_indices import *
from nrpy.infrastructures.ETLegacy.interface_ccl import construct_interface_ccl
from nrpy.helpers.colorize_text import colorize
from here import here

def flat_metric(out:Indexed, ni:Idx, nj:Idx)->Expr:
    i = to_num(ni)
    j = to_num(nj)
    if i==j:
        return sympify(1)
    else:
        return sympify(0)

# Create a set of grid functions
gf = GF()

# Declare gfs
p   = gf.decl("p",[li])
p_t = gf.decl("p_t",[li])
p_d = gf.decl("p_d",[li,lj])
u   = gf.decl("u",[])
u_t = gf.decl("u_t",[])
u_d = gf.decl("u_d",[ui])

siter2 = gf.decl("siter2", [li,lj])
gf.add_sym(siter2[li,lj], li, lj)
iter1 = gf.decl("iter1", [li])

# Declare the metric
g = gf.decl("g",  [li,lj])

# Declare params
spd = gf.add_param("spd", default=1.0, desc="The wave speed")

# Fill in values
gf.fill_in(g[li,lj], flat_metric)
here()
gf.fill_in(g[ui,uj], flat_metric)

# Fill in with defaults
here()
gf.fill_in(p[li], lambda _,i: mkSymbol(f"pD{to_num(i)}"))
here()
gf.fill_in(p_t[li], lambda _,i: mkSymbol(f"p_tD{to_num(i)}"))


# Fill in the deriv variables with a function call
#
gf.declfun("div1", True)
here()
gf.fill_in(p_d[li,lj], lambda _,i,j: div1(p[j],i))
here()
gf.fill_in(u_d[li], lambda _,i: div1(u,i))
here()
gf.fill_in(siter2[li,lj], alt=div1(p[li],lj), f=lambda _,i,j: mkSymbol(f"pD{to_num(i)}_dD{to_num(j)}"))
here()
gf.fill_in(iter1[lj], alt=div1(u,lj), f=lambda _,j: mkSymbol(f"u_dD{to_num(j)}"))

here()
res = gf.do_subs(spd*g[ui,uj]*div1(p[lj],li))

# Add the equations we want to evolve.
gf.add_eqn(p_t[lj], spd*div1(u,lj), "EVO")
gf.add_eqn(u_t, spd*g[ui,uj]*div1(p[lj],li), "EVO")

# Ensure the equations make sense
gf.diagnose()
exit(0)

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

#generate_Carpet_thorn(gf, "wavetoy", "WaveToyProj")
#generate_CarpetX_thorn(gf, "wavetoy", "WaveToyProj")
#generate_StandaloneRust(gf)
