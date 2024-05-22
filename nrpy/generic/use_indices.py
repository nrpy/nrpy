"""
Use the Sympy Indexed type for relativity expressions.
"""
from typing import Union, Set, Dict, List, Any, cast, Callable, Tuple, Optional, Type
from sympy import IndexedBase, Idx, Eq, Indexed, Basic, Mul, Expr, Eq, Symbol, Integer
from sympy.core.function import UndefinedFunction as UFunc
from inspect import currentframe
from nrpy.generic.sympywrap import *
from nrpy.generic.eqnlist import EqnList
from nrpy.generic.symm import Sym
from nrpy.helpers.colorize_text import colorize
import re
import sys

lookup_pair = dict()
def mkPair(s:str)->Tuple[Idx,Idx]:
    assert len(s)
    u, l = mkIdxs(f"u{s} l{s}")
    lookup_pair[l] = u
    lookup_pair[u] = l
    return u, l

def to_num(ind:Idx)->int:
    s = str(ind)
    assert s[0] in ["u", "l"]
    return int(s[1])

# Some basic indexes to use
ui, li = mkPair('i')
uj, lj = mkPair('j')
uk, lk = mkPair('k')
ua, la = mkPair('a')
ub, lb = mkPair('b')
uc, lc = mkPair('c')
u0, l0 = mkPair('0')
u1, l1 = mkPair('1')
u2, l2 = mkPair('2')
u3, l3 = mkPair('3')
u4, l4 = mkPair('4')
u5, l5 = mkPair('5')
up_indices = u0, u1, u2, u3, u4, u5 
down_indices = l0, l1, l2, l3, l4, l5 

multype = Mul #type(i*j)
addtype = type(ui + uj)
eqtype = Eq
powtype = type(ui**uj)

dimension = 3
def set_dimension(dim:int)->None:
    global dimension
    dimension = dim

ord0 = ord('0')
ord9 = ord('9')
def is_letter_index(sym:Basic)->bool:
    if type(sym) != Idx:
        return False
    s = str(sym)
    if len(s) != 2:
        return False
    if s[0] not in ["u","l"]:
        return False
    n = ord(s[1])
    return n < ord0 or n > ord9

def get_indices(xpr:Expr)->Set[Idx]:
    """ Return all indices of IndexedBase objects in xpr. """
    ret = set()
    for sym in finder(xpr):
        if is_letter_index(sym):
            ret.add(cast(Idx, sym))
    return ret
    ###
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

def get_free_indices(xpr:Expr)->Set[Idx]:
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

def get_contracted_indices(xpr:Expr)->Set[Idx]:
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
            sym_result = sym.apply(result)
            if result != sym_result:
                continue
        output += [(do_subs(xpr,index_values, sym),index_values.copy())]
    return output

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
    Callable[[Expr,Idx], Expr],
    Callable[[Expr,Idx,Idx], Expr],
    Callable[[Expr,Idx,Idx,Idx], Expr],
    Callable[[Expr,Idx,Idx,Idx,Idx], Expr],
    Callable[[Expr,Idx,Idx,Idx,Idx,Idx], Expr],
    Callable[[Expr,Idx,Idx,Idx,Idx,Idx,Idx], Expr]]

def fill_in_default_(out: Indexed, *inds:int)->Expr:
    return mksymbol_for_tensor(out)

fill_in_default = cast(fill_in_type, fill_in_default_)
            
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
        self.defn : Dict[str,Tuple[str,List[Idx]]] = dict()

    def add_param(self, name:str, default:param_default_type, desc:str, values:param_values_type=None)->Symbol:
        self.params[name] = Param(name, default, desc, values)
        return mkSymbol(name)

    def _add_eqn2(self, lhs2:Symbol, rhs2:Expr)->None:
        rhs2 = self.do_subs(expand_contracted_indices(rhs2, self.symmetries))
        if str(lhs2) in self.gfs:
            self.eqnlist.add_output(lhs2)
        for item in finder(rhs2):
            if str(item) in self.gfs:
                #assert item.is_Symbol
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
                lhs2 = cast(Symbol, self.do_subs(lhsx, self.subs))
                rhs2 = self.do_subs(rhs, inds, self.subs)
                rhs2 = self.do_subs(rhs2, inds, self.subs)
                self._add_eqn2(lhs2, rhs2)
        elif type(lhs) in [IndexedBase, Symbol]:
            lhs2 = cast(Symbol, self.do_subs(lhs, self.subs))
            eci = expand_contracted_indices(rhs, self.symmetries)
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
        self.symmetries.add(tens.base, i1, i2, sgn)
    
    def declfun(self, funname:str, is_stencil:bool)->UFunc:
        fun = mkFunction(funname)
        self.eqnlist.add_func(fun, is_stencil)

        # If possible, insert the symbol into the current environment
        frame = currentframe()
        f_back = None if frame is None else frame.f_back
        globs  = None if f_back is None else f_back.f_globals
        if globs is not None:
            globs[funname] = fun

        return fun

    def declscalar(self, basename:str)->Symbol:
        ret = mkSymbol(basename)
        self.gfs[basename] = ret
        self.defn[basename] = (basename, list())

        # If possible, insert the symbol into the current environment
        frame = currentframe()
        f_back = None if frame is None else frame.f_back
        globs  = None if f_back is None else f_back.f_globals
        if globs is not None:
            globs[basename] = ret

        return ret

    def decl(self, basename:str, indices:List[Idx])->IndexedBase:
        ret = mkIndexedBase(basename, shape=tuple([dimension]*len(indices)) )
        self.gfs[basename] = ret
        self.defn[basename] = (basename, list(indices))

        # If possible, insert the symbol into the current environment
        frame = currentframe()
        f_back = None if frame is None else frame.f_back
        globs  = None if f_back is None else f_back.f_globals
        if globs is not None:
            globs[basename] = ret

        return ret

    def show_tensortypes(self)->None:
        keys : Set[str] = set()
        for k1 in self.eqnlist.inputs:
            keys.add(str(k1))
        for k2 in self.eqnlist.outputs:
            keys.add(str(k2))
        for k in keys:
            group, indices, members = self.get_tensortype(k)
            print(colorize(k,"green"),"is a member of",colorize(group,"green"),"with indices",colorize(indices,"cyan"),"and members",colorize(members,"magenta"))

    def get_tensortype(self, item:Union[str,Math])->Tuple[str,List[Idx],List[str]]:
        k = str(item)
        assert k in self.gfs.keys(), f"Not a defined symbol {item}"
        v = self.base_of.get(k, None)
        if v is None:
            return("none", list(), list()) #scalar
        return (v, self.defn[v][1], self.groups[v])

    def fill_in(self, indexed:IndexedBase, f:fill_in_type=fill_in_default, alt:Any=None)->None:
        for tup in expand_free_indices(indexed, self.symmetries):
            out, indrep = tup
            assert type(out) == Indexed
            inds = out.indices
            subj : Expr
            if alt is None:
                subj = out
            else:
                subj = self.do_subs(alt, indrep)
            subval_ = f(subj, *inds)
                
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
            print(colorize(subj,"red"),colorize("->","magenta"),colorize(subval_,"cyan"))
            self.subs[subj] = subval_

    def expand_eqn(self, eqn:Eq)->List[Eq]:
        result : List[Eq] = list()
        for tup in expand_free_indices(eqn.lhs, sym):
            lhs, inds = tup
            result += [mkEq(self.do_subs(lhs, self.subs), self.do_subs(eqn.rhs, inds, self.subs))]
        return result

    #def expand(self, arg:Symbol)->Expr:
    #    return self.do_subs(expand_contracted_indices(arg, self.symmetries), self.subs)

    def do_subs(self, arg:Expr, *subs:do_subs_table_type)->Expr:
        for i in range(20):
            new_arg = arg
            new_arg = expand_contracted_indices(new_arg, self.symmetries)
            new_arg = do_subs(new_arg, self.subs, *subs)
            if new_arg == arg:
                return arg
            arg = new_arg
        raise Exception(str(arg))

if __name__ == "__main__":
    gf = GF()
    B = gf.decl("B",[lc,lb])
    M = gf.decl("M",[la,lb])

    # Anti-Symmetric
    gf.add_sym(M[la,lb], la, lb, -1)

    n = 0
    for out in gf.expand_eqn(mkEq(M[la,lb], B[la,lb])):
        print(out)
        n += 1
    assert n==3

    # Symmetric
    N = gf.decl("N",[la,lb])
    gf.add_sym(N[la,lb], la, lb, 1)

    n = 0
    for out in gf.expand_eqn(mkEq(N[la,lb], B[la,lb])):
        print(out)
        n += 1
    assert n==6

    # Non-Symmetric
    Q = gf.decl("Q",[la,lb])

    n = 0
    for out in gf.expand_eqn(mkEq(Q[la,lb], B[la,lb])):
        print(out)
        n += 1
    assert n==9
