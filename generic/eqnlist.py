from typing import TypeVar, Literal, List, Dict, Union, Tuple, Any, Set, Generic, Iterator, Optional, Type, cast, Any

from sympy.core.symbol import Symbol
from sympy.core.expr import Expr
from sympy import symbols, Function, diff
from sympy import cse as cse_ 

cse_return = Tuple[List[Tuple[Symbol,Expr]],List[Expr]]
def cse(arg:List[Expr])->cse_return:
    return cast(cse_return,cse_(arg)) # type: ignore[no-untyped-call]

class EqnList:
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

    def optimize(self)->None:
        self.trim()
        needed:Set[Symbol] = set()
        used:Set[Symbol] = set()
        temps:Set[Symbol] = set()
        self.order = list()
        for k in self.eqns:
            if k not in self.inputs and k not in self.params and k not in self.outputs:
                temps.add(k)
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
el.optimize()
el.dump()
el.cse()
el.optimize()
el.dump()
