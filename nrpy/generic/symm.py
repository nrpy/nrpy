from typing import Union, Set, Dict, List, Any, cast, Callable, Tuple, Optional, Type
from sympy import IndexedBase, Idx, Eq, Indexed, Basic, Mul, Expr, Eq, Symbol, Integer
from inspect import currentframe
from nrpy.generic.sympywrap import *
from nrpy.helpers.colorize_text import colorize

class Sym(Applier):
    def __init__(self)->None:
        self.sd : Dict[IndexedBase,List[Tuple[int,int,int]]] = dict()
        self.modified = False

    def add(self, m:IndexedBase, p1:int, p2:int, sgn:int)->None:
        assert p1 < p2
        assert sgn in [1, -1]
        if m not in self.sd:
            self.sd[m] = []
        self.sd[m].append((p1, p2, sgn))

    def match(self, expr:Basic)->bool:
        if not hasattr(expr, "base"):
            return False
        syms = self.sd.get(expr.base, None)
        if syms is None:
            return False
        for ind1, ind2, sgn in syms:
            s1 = str(expr.args[ind1+1])
            assert len(s1) == 2, f"s1='{s1}'"
            s2 = str(expr.args[ind2+1])
            assert len(s2) == 2, f"s2='{s2}'"
            if s1[0] != s2[0]:
                continue
            if s1 > s2:
                return True
            elif s1 == s2 and sgn < 0:
                return True

        return False

    def replace(self, expr:Indexed)->Any:
        syms = self.sd.get(expr.base, list())
        args = [cast(Idx, a) for a in expr.args[1:]]
        retsgn = 1
        for ind1, ind2, sgn in syms:
            s1 = str(args[ind1])
            s2 = str(args[ind2])
            if s1[0] != s2[0]:
                continue
            if s1 > s2:
                self.modified = True
                args[ind1], args[ind2] = args[ind2], args[ind1]
                retsgn *= sgn
            elif s1 == s2 and sgn < 0:
                self.modified = True
                return sympify(0)
        if retsgn == 1:
            ret =  expr.base.__getitem__(tuple(args))
        else:
            # It seems wrong that this is considered an Any
            ret = -expr.base.__getitem__(tuple(args))
        return ret

    def apply(self, expr:Basic)->Basic:
        self.modified = True
        while self.modified:
            self.modified = False
            expr = expr.replace(self.match, self.replace) # type: ignore[no-untyped-call]
        return expr

def test()->None:
    u0, u1, u2, u3 = mkIdxs('u0 u1 u2 u3')
    sym = Sym()
    eps = mkIndexedBase("eps",shape=(3,3,3))
    sym.add(eps,0,1,-1)
    sym.add(eps,0,2,-1)
    sym.add(eps,1,2,-1)
    assert sym.apply(eps[u2,u1,u0]) == -eps[u0,u1,u2]
    assert sym.apply(eps[u0,u2,u1]) == -eps[u0,u1,u2]
    assert sym.apply(eps[u0,u1,u2]) ==  eps[u0,u1,u2]
    assert sym.apply(eps[u0,u1,u1]) ==  0
if __name__ == "__main__":
    test()
