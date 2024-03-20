from typing import TypeVar, Literal, List, Dict, Union, Tuple, Any, Set, Generic, Iterator, Optional, Type, cast, Any

from sympy.core.symbol import Symbol
from sympy.core.expr import Expr
from sympy import symbols, Function, diff
from sympy import cse as cse_ 
from here import here

# Quieting mypy. Not sure if this trick is needed.
cse_return = Tuple[List[Tuple[Symbol,Expr]],List[Expr]]
def cse(arg:List[Expr])->cse_return:
    return cast(cse_return,cse_(arg)) # type: ignore[no-untyped-call]

class EqnList:
    """
    This class models a generic list of equations. As such, it knows nothing about the rest of NRPy+.
    Ultimately, the information in this class will be used to generate a loop to be output by NRPy+.
    All it knows are the following things:
    (1) params - These are quantities that are generated outside the loop.
    (2) inputs - These are quantities which are read by equations but never written by them.
    (3) outputs - These are quantities which are written by equations but never read by them.
    (4) equations - These relate inputs to outputs. These may contain temporary variables, i.e.
                    quantities that are both read and written by equations.

    This class can remove equations and parameters that are not needed, but will complain
    about inputs that are not needed. It can also detect errors in the classification of
    symbols as inputs/outputs/params.
    """
    def __init__(self)->None:
        self.eqns:Dict[Symbol,Expr] = dict()
        self.params:Set[Symbol] = set()
        self.inputs:Set[Symbol] = set()
        self.outputs:Set[Symbol] = set()
        self.order:List[Symbol] = list()

    def add_param(self, lhs:Symbol)->None:
        assert lhs not in self.outputs, f"The symbol '{lhs}' is alredy in outputs"
        assert lhs not in self.inputs, f"The symbol '{lhs}' is alredy in outputs"
        self.params.add(lhs)

    def add_input(self, lhs:Symbol)->None:
        assert lhs not in self.outputs, f"The symbol '{lhs}' is alredy in outputs"
        assert lhs not in self.params, f"The symbol '{lhs}' is alredy in outputs"
        self.inputs.add(lhs)

    def add_output(self, lhs:Symbol)->None:
        assert lhs not in self.inputs, f"The symbol '{lhs}' is alredy in outputs"
        assert lhs not in self.params, f"The symbol '{lhs}' is alredy in outputs"
        self.outputs.add(lhs)

    def add_eqn(self, lhs:Symbol, rhs:Expr)->None:
        self.eqns[lhs] = rhs

    def diagnose(self)->None:
        """ Discover inconsistencies and errors in the param/input/output/equation sets. """
        needed:Set[Symbol] = set()
        used:Set[Symbol] = set()
        temps:Set[Symbol] = set()
        self.order = list()

        read:Set[Symbol] = set()
        written:Set[Symbol] = set()

        for k in self.eqns:
            written.add(k)
            for q in self.eqns[k].free_symbols:
                read.add(q)

        for k in written:
            if k in read:
                temps.add(k)
            else:
                assert k in self.outputs, f"Symbol '{k}' is written, but it is not a temp and not in outputs"

        for k in read:
            assert k in self.inputs or self.params or temps, f"Symbol '{k}' is read, but it is not a temp, parameter, or input."

        for k in self.outputs:
            needed.add(k)
            assert k in self.eqns, f"The symbol '{k}' appears in outputs, but no eqn assigns to it."
        again = True
        wait = False
        while again:
            again = False
            for k,v in enumerate(self.eqns):
                assert k not in self.inputs, f"The symbol '{k}' appears in inputs, but it is assigned to by '{k} = {v}'"
                assert k not in self.params, f"The symbol '{k}' appears in params, but it is assigned to by '{k} = {v}'"
                if k in needed and k not in used:
                    if k in temps or k in self.outputs:
                        wait = False
                        for k2 in v.free_symbols:
                            if k2 in temps and k2 not in used:
                                wait = True
                                if k2 not in needed:
                                    needed.add(cast(Symbol,k2))
                                    again = True
                    if wait:
                        continue
                    used.add(k)
                    self.order.append(k)
                    for k2 in v.free_symbols:
                        k3 = cast(Symbol, k2)
                        needed.add(k3)
                        assert k3 not in self.outputs, f"The symbol '{k3}' appears in outputs, but it is read by '{k} = {v}'"
                    again = True
        assert not wait, "Cycle in temp assignment"
        for k in needed:
            assert k in self.inputs or k in self.params or k in self.eqns,\
                f"The symbol '{k}' is needed but is not defined"
        for k in self.inputs:
            if k not in needed:
                print(f"Warning: symbol '{k}' appears in inputs but is not needed")

    def trim(self)->None:
        """ Remove temporaries of the form "a=b". They are clutter. """
        subs:Dict[Symbol,Symbol] = dict()
        for k,v in self.eqns.items():
            if v.is_symbol:
                # k is not not needed
                subs[k] = v
                print(f"Warning: equation '{k} = {v}' can be trivially eliminated")

        new_eqns:Dict[Symbol,Expr] = dict()
        for k in self.eqns:
            if k not in subs:
                v = self.eqns[k]
                v2 = v.subs(subs)
                new_eqns[k] = v2

        self.eqns = new_eqns

    def cse(self)->None:
        """ Invoke Sympy's CSE method, but ensure that the order of the resulting assignments is correct. """
        indexes:List[Symbol]=list()
        old_eqns:List[Expr]=list()
        for k in self.eqns:
            indexes.append(k)
            old_eqns.append(self.eqns[k])
        new_eqns, mod_eqns = cse(old_eqns) 
        e : Tuple[Symbol, Expr]
        for e in new_eqns:
            assert e[0] not in self.inputs and e[0] not in self.params and e[0] not in self.eqns
            self.add_eqn(e[0], e[1])
        for i in range(len(indexes)):
            k = indexes[i]
            v = old_eqns[i]
            m = mod_eqns[i]
            self.eqns[k] = m

    def dump(self)->None:
        print("Dumping Equations:")
        for k in self.order:
            print(" ",k,"->",self.eqns[k])

if __name__ == "__main__":
    a, b, c, d, e, f, g = symbols("a b c d e f g")
    el = EqnList()
    el.add_input(a)
    el.add_input(f)
    el.add_input(b)
    el.add_output(d)
    el.add_eqn(c, e+b+a**3)
    el.add_eqn(c, g+f+a**3)
    el.add_eqn(e, a**2)
    el.add_eqn(g, e)
    el.add_eqn(d, c*b+a**3)
    el.diagnose()
    el.dump()
    el.cse()
    el.diagnose()
    el.dump()