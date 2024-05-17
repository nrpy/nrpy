from typing import Tuple, List, Dict, Any, Union, cast
from sympy import cse as cse_, IndexedBase, Idx, Symbol, Expr, Eq, Basic, sympify as sympify_, Expr, Mul, Indexed, Function
import re
from abc import ABC, abstractmethod
from sympy.core.function import UndefinedFunction as UFunc

class Applier(ABC):

    @abstractmethod
    def apply(self, Basic:Basic)->Basic:
        ...

IndexType = Union[Idx,Mul]

cse_return = Tuple[List[Tuple[Symbol,Expr]],List[Expr]]
def cse(arg:List[Expr])->cse_return:
    return cast(cse_return, cse_(arg)) # type: ignore[no-untyped-call]

def mkIdx(name:str)->Idx:
    return Idx(name) # type: ignore[no-untyped-call]

def mkSymbol(name:str)->Symbol:
    return Symbol(name) # type: ignore[no-untyped-call]

def mkFunction(name:str)->UFunc:
    return Function(name) # type: ignore[no-any-return]

def mkEq(a:Basic, b:Basic)->Eq:
    return Eq(a,b) # type: ignore[no-untyped-call]

def mkIdxs(names:str)->Tuple[Idx,...]:
    return tuple([Idx(name) for name in re.split(r'\s+', names)]) # type: ignore[no-untyped-call]

def mkIndexedBase(basename:str, shape:Tuple[int,...])->IndexedBase:
    return IndexedBase(basename, shape=shape) # type: ignore[no-untyped-call]

def mkIndexed(base:IndexedBase, *args:Union[int,IndexType])->Indexed:
    return Indexed(base, *args) # type: ignore[no-untyped-call]

def sympify(arg:Any)->Expr:
    return cast(Expr, sympify_(arg)) # type: ignore[no-untyped-call]

do_subs_table_type = Union[
        Dict[Idx,Idx],
        Dict[Indexed, Indexed],
        Dict[Expr, Expr],
        Dict[Symbol, Symbol],
        Applier
        ]

def do_subs(sym:Expr, *tables:do_subs_table_type)->Expr:
    result = sym
    for table in tables:
        if isinstance(table, Applier):
            result = cast(Expr, table.apply(result))
        else:
            result = cast(Expr, result.subs(table)) # type: ignore[no-untyped-call]
    return result
