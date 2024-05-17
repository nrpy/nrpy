"""
Use the Sympy Indexed type for relativity expressions.
"""
from typing import Union, Set, Dict, List, Any, cast, Callable, Tuple, Optional, Type
from sympy import IndexedBase, Idx, Eq, Indexed, Basic, Mul, Expr, Eq, Symbol, Integer
from inspect import currentframe
from nrpy.generic.sympywrap import *
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.symm import Sym
from nrpy.helpers.colorize_text import colorize
import re
import sys
from here import here

lookup_pair = dict()
def mkPair(s:str)->Tuple[Idx,Idx]:
    assert len(s)
    u, l = mkIdxs(f"u{s} l{s}")
    lookup_pair[l] = u
    lookup_pair[u] = l
    return u, l

# Some basic indexes to use
ui, li = mkPair('i')
uj, lj = mkPair('j')
uk, lk = mkPair('k')
ua, la = mkPair('a')
ub, lb = mkPair('b')
uc, lc = mkPair('c')
up_indices = mkIdxs('u0 u1 u2 u3 u4 u5')
down_indices = mkIdxs('l0 l1 l2 l3 l4 l5')
u0, u1, u2, u3, u4, u5 = up_indices
l0, l1, l2, l3, l4, l5 = down_indices

multype = Mul #type(i*j)
addtype = type(ui + uj)
eqtype = Eq
powtype = type(ui**uj)

dimension = 3
def set_dimension(dim:int)->None:
    global dimension
    dimension = dim

def get_indices(xpr:Basic)->Set[Idx]:
    """ Return all indices of IndexedBase objects in xpr. """
    ret = set()
    if type(xpr) in [multype, addtype, powtype]:
        for arg in xpr.args:
            ret.update(get_indices(arg))
        return ret
    elif hasattr(xpr, "indices"):
        ret.update(xpr.indices)
    return ret

def byname(x:Idx)->str:
    """ Return a string suitable for sorting a list of upper/lower indices. """
    s = str(x)
    assert len(s) == 2
    return s[1] + s[0]

num0 = ord('0')
num9 = ord('9')
def is_numeric_index(x:Idx)->bool:
    s = str(x)
    assert len(s) == 2
    n = ord(s[1])
    return num0 <= n and n <= num9

def get_numeric_index_value(x:Idx)->int:
    s = str(x)
    assert is_numeric_index(x)
    return ord(s[1]) - num0

def is_lower(x:Idx)->bool:
    s = str(x)
    return s[0] == 'l'

def is_upper(x:Idx)->bool:
    s = str(x)
    return s[0] == 'u'

def get_pair(x:Idx)->Tuple[Idx, Idx]:
    if is_lower(x):
        return lookup_pair[x], x
    else:
        return x, lookup_pair[x]

def is_pair(a:Idx, b:Idx)->bool:
    sa = str(a)
    sb = str(b)
    assert len(sa) == 2
    assert len(sb) == 2
    if sa[1] == sb[1] and ((sa[0] == 'u' and sb[0] == 'l') or (sa[0] == 'l' and sb[0] == 'u')):
        return True
    else:
        return False

# Check that this works
assert is_pair(ui, li)
assert is_pair(li, ui)
assert not is_pair(ui, lj)
assert not is_pair(li, uj)

def get_free_indices(xpr:Basic)->Set[Idx]:
    """ Return all uncontracted indices in xpr. """
    indices = list(get_indices(xpr))
    indices = sorted(indices, key=byname)
    ret = set()
    i = 0
    while i < len(indices):
        if i+1 < len(indices) and is_pair(indices[i], indices[i+1]):
            i += 2
        else:
            ret.add(indices[i])
            i += 1
    return ret

M = mkIndexedBase('M',(3,3))
assert sorted(list(get_free_indices(M[ui,uj]*M[lj, lk])), key=byname) == [ui, lk]

def get_contracted_indices(xpr:Basic)->Set[Idx]:
    """ Return all contracted indices in xpr. """
    indices = list(get_indices(xpr))
    indices = sorted(indices, key=byname)
    ret = set()
    i = 0
    while i < len(indices):
        if i+1 < len(indices) and is_pair(indices[i], indices[i+1]):
            ret.add(indices[i])
            i += 2
        else:
            i += 1
    return ret

assert sorted(list(get_contracted_indices(M[ui,uj]*M[lj, lk])), key=byname) == [lj]

def incr(index_list:List[Idx], index_values:Dict[Idx, Idx])->bool:
    """ Increment the indices in index_list, creating an index_values table with all possible permutations. """
    if len(index_list) == 0:
        return False
    ix = 0
    if len(index_values)==0:
        for ind_ in index_list:
            uind, ind = get_pair(ind_)
            index_values[ind] = l0
            index_values[uind] = u0
        return True
    while True:
        if ix >= len(index_list):
            return False
        uind, ind = get_pair(index_list[ix])
        index_value = get_numeric_index_value(index_values[ind])
        if index_value == dimension-1:
            index_values[ind] = l0
            index_values[uind] = u0
            ix += 1
        else:
            index_values[ind] = down_indices[index_value+1]
            index_values[uind] = up_indices[index_value+1]
            break
    return True

# Check
ilist = [ui, lj]
dvals : Dict[Idx,Idx] = dict()
valscount = 0
while incr(ilist, dvals):
    valscount += 1
assert valscount == dimension**2

def expand_contracted_indices(xpr:Expr, sym:Sym)->Expr:
    if type(xpr) == addtype:
        ret : Expr = sympify(0)
        for arg in xpr.args:
            ret += expand_contracted_indices(arg, sym)
        return ret
    index_list = sorted(list(get_contracted_indices(xpr)), key=byname)
    if len(index_list) == 0:
        return xpr
    output = sympify(0)
    index_values : Dict[Idx, Idx]= dict()
    while incr(index_list, index_values):
        output += do_subs(xpr, index_values, sym)
    return output

# Check
sym = Sym()
assert expand_contracted_indices(M[ui, li], sym) == M[u0, l0] + M[u1, l1] + M[u2, l2]
assert expand_contracted_indices(M[ui, lj]*M[li, uk], sym) == M[l0, uk]*M[u0, lj] + M[l1, uk]*M[u1, lj] + M[l2, uk]*M[u2, lj]

def expand_free_indices(xpr:Expr, sym:Sym)->List[Tuple[Expr, Dict[Idx, Idx]]]:
    index_list = sorted(list(get_free_indices(xpr)), key=str)
    output : List[Tuple[Expr, Dict[Idx, Idx]]] = list()
    xpr = expand_contracted_indices(xpr, sym)
    index_values : Dict[Idx, Idx] = dict()
    while incr(index_list, index_values):
        assert len(index_values) != 0, "Something very bad happened"
        if type(xpr) == Indexed:
            result = do_subs(xpr, index_values)
            if result == sym.apply(result):
                continue
        output += [(do_subs(xpr,index_values, sym),index_values.copy())]
    return output

## def add_asym(tens:Indexed, ix1:Idx, ix2:Idx)->None:
##     """
##     add_asym(tens[ua,ub,uc,ud], uc, ud) would make a tensor
##     that is antisymmetric in the last two indices.
##     """
##     add_sym(tens, ix1, ix2, sgn=-1)

## def add_sym(tens:Indexed, ix1:Idx, ix2:Idx, sgn:int=1)->None:
##     """
##     add_sym(tens[ua,ub,uc,ud], uc, ud) would make a tensor
##     that is symmetric in the last two indices.
##     """
##     assert type(tens) == Indexed
##     assert type(tens.args[0]) == IndexedBase, f"tens.args[0]={type(tens.args[0])}"
##     base:IndexedBase = tens.args[0]
##     i1 = -1
##     i2 = -1
##     for i in range(1,len(tens.args)):
##         if tens.args[i] == ix1:
##             i1 = i-1
##         if tens.args[i] == ix2:
##             i2 = i-1
##     assert i1 != -1, f"Index {ix1} not in {tens}"
##     assert i2 != -2, f"Index {ix2} not in {tens}"
##     assert i1 != i2, f"Index {ix1} cannot be symmetric with itself in {tens}"
##     if i1 > i2:
##         i1, i2 = i2, i1
##     index_list : List[Idx] = cast(List[Idx], list(tens.args)[1:])
##     # create an index list with i1 and i2 swapped
##     index_list2 = \
##         index_list[:i1] + \
##         index_list[i2:i2+1] + \
##         index_list[i1+1:i2] + \
##         index_list[i1:i1+1] + \
##         index_list[i2+1:]
##     index_values : Dict[Idx,Idx] = dict()
##     while incr(index_list, index_values):
##         if get_numeric_index_value(index_values[ix1]) > get_numeric_index_value(index_values[ix2]):
##             args1 = [index_values[ix] for ix in index_list]
##             args2 = [index_values[ix] for ix in index_list2]
##             term1 = mkIndexed(base, *args1)
##             term2 = mkIndexed(base, *args2)
##             symmetries[term1] = sgn*term2
##             symmetries[term1] = sgn*term2

## # Check
## add_sym(M[da,db],da,db)
## assert M[d0,d1].subs(symmetries) == M[d1,d0].subs(symmetries)
## # It would be bad if anyone actually did this...
## # A tensor can't really be antisymmetric in its upper indices
## # and symmetric in its lower indices.
## add_asym(M[ua,ub],ua,ub)
## assert M[u0,u1].subs(symmetries) == -M[u1,u0].subs(symmetries)

def mksymbol_for_tensor(out:Indexed)->Symbol:
    """
    Define a symbol for a tensor using standard NRPy+ rules.
    For an upper index put a U, for a lower index put a D.
    Follow the string of U's and D's with the integer value
    of the up/down index.

    :param out: The tensor expression with integer indices.

    :return: a new sympy symbol
    """
    base:str = str(out.args[0])
    for i in out.args[1:]:
        assert type(i) == Idx
        if is_upper(i):
            base += "U"
        else:
            base += "D"
    for i in out.args[1:]:
        assert type(i) == Idx
        base += str(get_numeric_index_value(i))
    return mkSymbol(base)

# It's horrible that Python can't let me type this any other way
fill_in_type = Union[
    Callable[[Indexed,int], Expr],
    Callable[[Indexed,int,int], Expr],
    Callable[[Indexed,int,int,int], Expr],
    Callable[[Indexed,int,int,int,int], Expr],
    Callable[[Indexed,int,int,int,int,int], Expr],
    Callable[[Indexed,int,int,int,int,int,int], Expr]]

def fill_in_default_(out: Indexed, *inds:int)->Expr:
    return mksymbol_for_tensor(out)

fill_in_default = cast(fill_in_type, fill_in_default_)

## subs : Dict[Expr, Expr] = dict()
## def fill_in(indexed:IndexedBase, f:fill_in_type=fill_in_default)->None:
##     for tup in expand_free_indices(indexed):
##         out, _ = tup
##         assert type(out) == Indexed
##         inds = out.indices
##         subs[out] = f(out, *inds)

## # Check
## fill_in(M[da,db])
## assert len(subs) == 6 # M is symmetric
## subs = dict()

## Q = mkIndexedBase('Q',(3,3))
## fill_in(Q[da,db])
## assert len(subs) == 9 # Q has no symmetries
## subs = dict()

#def expand(arg:Expr, sym:Sym)->Expr:
#    return do_subs(expand_contracted_indices(arg, sym), subs)
            
param_default_type = Union[float,int,str,bool]
param_values_type = Optional[Union[Tuple[float,float],Tuple[int,int],Tuple[bool,bool],str,Set[str]]]
min_max_type = Union[Tuple[float,float],Tuple[int,int]]

class Param:
    def __init__(self, name:str, default:param_default_type, desc:str, values:param_values_type)->None:
        self.name = name
        self.values = values
        self.desc = desc
        self.default = default

    def get_min_max(self)->min_max_type:
        ty = self.get_type()
        if ty == int:
            if self.values is not None:
                return cast(min_max_type, self.values)
            return (-2**31, 2**31-1)
        elif ty == float:
            if self.values is not None:
                return cast(min_max_type, self.values)
            return (sys.float_info.min, sys.float_info.max)
        else:
            assert False

    def get_values(self)->param_values_type:
        if self.values is not None:
            return self.values
        ty = self.get_type()
        if ty == bool:
            return (False, True)
        elif ty == str:
            return ".*"
        else:
            return self.get_min_max()

    def get_type(self)->Type[Any]:
        if self.values == None:
            return type(self.default)
        elif type(self.values) == set:
            assert type(self.default) == str
            return set # keywords
        elif type(self.values) == str:
            # values is a regex
            assert type(self.default) == str
            return str
        elif type(self.values) == tuple and len(self.values) == 2:
            assert type(self.default) in [int, float]
            assert type(self.values[0]) in [int, float]
            assert type(self.values[1]) in [int, float]
            if type(self.default)==float or type(self.values[0]) == float or type(self.values[1]) == float:
                return float
            else:
                return int
        else:
            assert False

####
div1 = mkFunction("div1")
val = mkSymbol("val")
x = mkSymbol("x")
y = mkSymbol("y")
z = mkSymbol("z")

#def div1match(u:Symbol)->bool:
#    if hasattr(u, "func"):
#        fn = getattr(u, "func")
#        if fn == div1:
#            return True
#    return False

#def div1repl(fun:Symbol)->Expr:
#    here("fun:",fun, fun.args)
#    funstr = str(fun.args[0])
#    g = re.match(r"(.*_d[DU]+)(\d+)$", funstr)
#    if g:
#        funstr = g.group(1)
#        divargs = g.group(2)
#    else:
#        funstr += "_d"
#        divargs = ""
#    arg2 = fun.args[1]
#    if arg2 == x:
#        arg2 = -1
#    elif arg2 == y:
#        arg2 = -2
#    elif arg2 == z:
#        arg2 = -3
#    else:
#        here("arg2:",arg2)
#        arg2 = int(arg2)
#    assert arg2 != 0, "For div1repl, we want the base index to be 1 not zero"
#    if arg2 > 0:
#        funstr += "U"
#    else:
#        funstr += "D"
#    divargs += str(abs(arg2) - 1)
#    funstr += divargs
#    return Symbol(funstr)

####

class GF:
    def __init__(self)->None:
        self.symmetries = Sym()
        self.gfs:Dict[str,Union[Indexed, IndexedBase, Symbol]] = dict()
        self.subs : Dict[Expr, Expr] = dict()
        self.eqnlist : EqnList = EqnList()
        self.params : Dict[str,Param] = dict()
        self.base_of : Dict[str, str] = dict()
        self.groups : Dict[str, List[str]] = dict()
        self.props : Dict[str,List[Integer]] = dict()
        self.defn : Dict[str,str] = dict()

    def add_param(self, name:str, default:param_default_type, desc:str, values:param_values_type=None)->Symbol:
        self.params[name] = Param(name, default, desc, values)
        return mkSymbol(name)

    def _add_eqn2(self, lhs2:Symbol, rhs2:Expr)->None:
        if str(lhs2) in self.gfs:
            self.eqnlist.add_output(lhs2)
        for item in rhs2.free_symbols:
            if str(item) in self.gfs:
                assert item.is_Symbol
                self.eqnlist.add_input(cast(Symbol, item))
            elif str(item) in self.params:
                assert item.is_Symbol
                self.eqnlist.add_param(cast(Symbol, item))
        self.eqnlist.add_eqn(lhs2, rhs2)

    def add_eqn(self, lhs:Union[Indexed,IndexedBase,Symbol], rhs:Symbol, eqntype:str)->None:
        lhs2 : Symbol
        if type(lhs) == Indexed:
            for tup in expand_free_indices(lhs, self.symmetries):
                lhsx, inds = tup
                here("inds:",inds)
                here("lhsx:",lhsx)
                lhs2 = cast(Symbol, self.do_subs(lhsx, self.subs))
                here("lhsx:",lhsx)
                here("rhs:",rhs)
                rhs2 = self.do_subs(rhs, inds, self.subs)
                here("rhs:",rhs)
                self._add_eqn2(lhs2, rhs2)
        elif type(lhs) in [IndexedBase, Symbol]:
            lhs2 = cast(Symbol, self.do_subs(lhs, self.subs))
            eci = expand_contracted_indices(rhs, self.symmetries)
            here("eci:",eci, self.subs)
            rhs2 = self.do_subs(eci, self.subs)
            self._add_eqn2(lhs2, rhs2)
        else:
            print("other:",lhs,rhs,type(lhs),type(rhs))
            raise Exception()

    def cse(self)->None:
        self.eqnlist.cse()

    def dump(self)->None:
        self.eqnlist.dump()

    def diagnose(self)->None:
        self.eqnlist.diagnose()

    def add_sym(self, tens:Indexed, ix1:Idx, ix2:Idx, sgn:int=1)->None:
        assert type(tens) == Indexed
        base : IndexedBase = cast(IndexedBase, tens.args[0])
        i1 = -1
        i2 = -1
        for i in range(1,len(tens.args)):
            if tens.args[i] == ix1:
                i1 = i-1
            if tens.args[i] == ix2:
                i2 = i-1
        assert i1 != -1, f"Index {ix1} not in {tens}"
        assert i2 != -2, f"Index {ix2} not in {tens}"
        assert i1 != i2, f"Index {ix1} cannot be symmetric with itself in {tens}"
        if i1 > i2:
            i1, i2 = i2, i1
        sym.add(tens.base, i1, i2, sgn)
    
    def decl(self, basename:str, indices:List[Idx])->IndexedBase:
        frame = currentframe()
        f_back = None if frame is None else frame.f_back
        globs  = None if f_back is None else f_back.f_globals
        ret = mkIndexedBase(basename, shape=tuple([dimension]*len(indices)) )
        self.gfs[basename] = ret
        self.defn[basename] = f"{basename}{indices}"
        if globs is not None:
            globs[basename] = ret
        return ret

    def fill_in(self, indexed:IndexedBase, f:fill_in_type=fill_in_default, base_zero:bool=True)->None:
        for tup in expand_free_indices(indexed, self.symmetries):
            out, _ = tup
            assert type(out) == Indexed
            inds = out.indices
            if base_zero:
                inds = [abs(i)-1 for i in inds]
            subval_ = f(out, *inds)
            if subval_.is_Number:
                pass
            elif subval_.is_Function:
                pass
            else:
                assert subval_.is_Symbol, f"{type(subval_)}, {subval_.__class__}, {subval_.is_Function}"
                subval = cast(Symbol, subval_)
                self.gfs[str(subval)] = subval
                self.base_of[str(subval)] = str(out.base)
                if str(out.base) not in self.groups:
                    self.groups[str(out.base)] = list()
                members = self.groups[str(out.base)]
                members.append(str(subval))
                print(f"base of {subval} is {out.base}")
                print(f"groups of {out.base} is {members}")
            print(colorize(subval_,"red"),colorize("->","magenta"),colorize(inds,"cyan"))
            self.subs[out] = subval_
            self.props[str(subval_)] = out.indices

    def expand_eqn(self, eqn:Eq)->List[Eq]:
        result : List[Eq] = list()
        for tup in expand_free_indices(eqn.lhs, sym):
            lhs, inds = tup
            result += [mkEq(self.do_subs(lhs, self.subs), self.do_subs(eqn.rhs, inds, self.subs))]
        return result

    def expand(self, arg:Symbol)->Expr:
        return self.do_subs(expand_contracted_indices(arg, self.symmetries), self.subs)

    def do_subs(self, arg:Expr, *subs:do_subs_table_type)->Expr:
        return do_subs(arg, *subs)

if __name__ == "__main__":
    gf = GF()
    gf.decl("M",[la,lb])
    gf.add_sym(M[la,lb], la, lb)
    gf.decl("B",[lc,lb])

    B : IndexedBase
    for out in gf.expand_eqn(mkEq(M[la,lb], B[la,lb])):
        print(out)
    pass
