"""
Use the Sympy Indexed type for relativity expressions.
"""
from sympy import symbols, IndexedBase, Idx, sympify, Eq, Indexed
from here import here
from inspect import currentframe

i, j, k = symbols('i j k', cls=Idx)
multype = type(i*j)
addtype = type(i+j)
eqtype = type(Eq(i,j))
powtype = type(i**j)

dimension = 3
def set_dimension(dim):
    global dimension
    dimension = dim

symmetries = dict()

def get_indices(xpr):
    """ Return all indices of IndexedBase objects in xpr. """
    ret = set()
    if type(xpr) in [multype, addtype, powtype]:
        for arg in xpr.args:
            ret.update(get_indices(arg))
        return ret
    elif hasattr(xpr, "indices"):
        ret.update(xpr.indices)
    return ret

def byname(x):
    """ Return a string suitable for sorting a list of upper/lower indices. Use negative indices as down indices. """
    s = str(x)
    if s[0] == "-":
        return s[1:]+"-"
    else:
        return s

def get_free_indices(xpr):
    """ Return all uncontracted indices in xpr. """
    indices = list(get_indices(xpr))
    indices = sorted(indices, key=byname)
    ret = set()
    i = 0
    while i < len(indices):
        if i+1 < len(indices) and indices[i] == -indices[i+1]:
            i += 2
        else:
            ret.add(indices[i])
            i += 1
    return ret

def get_contracted_indices(xpr):
    """ Return all contracted indices in xpr. """
    indices = list(get_indices(xpr))
    indices = sorted(indices, key=byname)
    ret = set()
    i = 0
    while i < len(indices):
        if i+1 < len(indices) and indices[i] == -indices[i+1]:
            ret.add(indices[i])
            i += 2
        else:
            i += 1
    return ret

def incr(index_list, index_values):
    """ Increment the indices in index_list, creating an index_values table with all possible permutations. """
    if len(index_list) == 0:
        return False
    ix = 0
    if len(index_values)==0:
        for ind_ in index_list:
            if type(ind_) != Idx:
                ind = -ind_
            else:
                ind = ind_
            index_values[ind] = 1
            index_values[-ind] = -1
        return True
    while True:
        if ix >= len(index_list):
            return False
        ind = index_list[ix]
        if type(ind) != Idx:
            ind = -ind
        if index_values[ind] == dimension:
            index_values[ind] = 1
            index_values[-ind] = -1 #index_values[ind]
            ix += 1
        else:
            index_values[ind] += 1
            index_values[-ind] = -index_values[ind]
            break
    return True

def expand_contracted_indices(xpr):
    if type(xpr) == addtype:
        ret = sympify(0)
        for arg in xpr.args:
            ret += expand_contracted_indices(arg)
        return ret
    index_list = sorted(list(get_contracted_indices(xpr)), key=str)
    if len(index_list) == 0:
        return xpr
    output = sympify(0)
    index_values = dict()
    while incr(index_list, index_values):
        output += xpr.subs(index_values).subs(symmetries)
    return output

def expand_free_indices(xpr):
    index_list = sorted(list(get_free_indices(xpr)), key=str)
    output = list()
    if type(xpr) == addtype:
        new_xpr = sympify(0)
        for arg in xpr.args:
            arg_index_list = sorted(list(get_free_indices(arg)), key=str)
            assert arg_index_list == index_list, f"Invalid tensor expression '{xpr}'"
            new_xpr += expand_contracted_indices(arg)
        xpr = new_xpr
    elif type(xpr) == multype or type(xpr) == Indexed:
        xpr = expand_contracted_indices(xpr)
    elif type(xpr) == eqtype:
        xpr = Eq(xpr.lhs, expand_contracted_indices(xpr.rhs))
        index_list = sorted(list(get_free_indices(xpr.lhs)), key=str)
    else:
        assert False, f"Type is {type(xpr)}"
    index_values = dict()
    while incr(index_list, index_values):
        assert len(index_values) != 0, "Something very bad happened"
        if type(xpr) == Indexed and xpr.subs(index_values) in symmetries:
            continue
        if type(xpr) == eqtype and xpr.lhs.subs(index_values) in symmetries:
            continue
        output += [xpr.subs(index_values).subs(symmetries)]
    return output

def add_asym(tens, ix1, ix2):
    return add_sym(tens, ix1, ix2, sgn=-1)

def add_sym(tens, ix1, ix2, sgn=1):
    assert type(tens) == Indexed
    base = tens.args[0]
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
    index_list = list(tens.args)[1:]
    # create an index list with i1 and i2 swapped
    index_list2 = \
        index_list[:i1] + \
        index_list[i2:i2+1] + \
        index_list[i1+1:i2] + \
        index_list[i1:i1+1] + \
        index_list[i2+1:]
    index_values = dict()
    while incr(index_list, index_values):
        if index_values[ix1] > index_values[ix2]:
            args1 = [index_values[ix] for ix in index_list]
            args2 = [index_values[ix] for ix in index_list2]
            term1 = Indexed(base, *args1)
            term2 = Indexed(base, *args2)
            symmetries[term1] = sgn*term2

def mksymbol_for_tensor(out):
    """
    Define a symbol for a tensor using standard NRPy+ rules.
    For an upper index put a U, for a lower index put a D.
    Follow the string of U's and D's with the integer value
    of the up/down index.

    :param out: The tensor expression with integer indices.

    :return: a new sympy symbol
    """
    base = str(out.args[0])
    for i in out.args[1:]:
        if i > sympify(0):
            base += "U"
        else:
            base += "D"
    for i in out.args[1:]:
        base += str(abs(i))
    return symbols(base)

subs = dict()
def fill_in(indexed:IndexedBase, f, normalize:bool=True):
    for out in expand_free_indices(indexed):
        inds = out.indices
        if normalize:
            inds = tuple([abs(i)-1 for i in list(inds)])
        subs[out] = f(*inds)

def expand(arg):
    return expand_contracted_indices(arg).subs(subs)

IndexType = Union[Idx,multype]
            
class GF:
    def __init__(self):
        self.symmetries = dict()

    def add_sym(self, tens:Indexed, ix1:IndexType, ix2:IndexType):
        assert type(tens) == Indexed
        base = tens.args[0]
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
        index_list = list(tens.args)[1:]
        # create an index list with i1 and i2 swapped
        index_list2 = \
            index_list[:i1] + \
            index_list[i2:i2+1] + \
            index_list[i1+1:i2] + \
            index_list[i1:i1+1] + \
            index_list[i2+1:]
        index_values = dict()
        while incr(index_list, index_values):
            if index_values[ix1] > index_values[ix2]:
                args1 = [index_values[ix] for ix in index_list]
                args2 = [index_values[ix] for ix in index_list2]
                term1 = Indexed(base, *args1)
                term2 = Indexed(base, *args2)
                self.symmetries[term1] = sgn*term2
    
    def decl(self, *args):
        globs = currentframe().f_back.f_globals
        for i in range(0,len(args),2):
            basename = args[i]
            indices = args[i+1]
            globs[basename] = IndexedBase(basename, shape=(dimension)*len(indices) )

if __name__ == "__main__":
    #===========================
    # IndexedBase is a tensor object without indices
    # Indexed is a tensor object with indices
    # indices are of type Idx
    M = IndexedBase('M')
    P = IndexedBase('P')
    R = IndexedBase('R')
    i, j, k = symbols('i j k', cls=Idx)

    # Negative index is contracted
    mm = M[i,j]*M[-j,k]

    # Define symmetries
    add_asym(P[i,j],i,j)
    add_asym(M[i,j],i,j)

    # Want to generate code for an equation like this:
    form = Eq(P[i,k], mm + R[i,k])

    print("form:",form.lhs,"=",form.rhs)

    # Create substitution rules for code generation
    mysubs = dict()
    for out in expand_free_indices(M[-i,j]) + \
            expand_free_indices(M[i,j]) + \
            expand_free_indices(R[i,j]) + \
            expand_free_indices(P[i,j]):
        # mksymbol_for_tensor uses the standard NRPy+ notation
        mysubs[out] = mksymbol_for_tensor(out)

    # Generate the code expressions by iterating over
    # the free (non-contracted) indices in the equation
    for out in  expand_free_indices(form):
        # Substitute the resluts with mysubs
        # Internally, the rhs will have contracted indices expanded
        print(out.lhs.subs(mysubs),"=",out.rhs.subs(mysubs))
