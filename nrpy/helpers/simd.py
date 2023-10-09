"""
Convert Expression to SIMD Intrinsics.

Authors: Ken Sible, Thiago Assumpcao, and Zachariah Etienne
Emails: ksible *at* outlook *dot** com
        assumpcaothiago *at* gmail *dot** com
        zachetie *at* gmail *dot** com
"""
from typing import Dict, Union, Any, Optional
import shutil
from pathlib import Path
from sympy import (
    Basic,
    Expr,
    Float,
    Function,
    Integer,
    Rational,
    Symbol,
    Abs,
    Add,
    cos,
    exp,
    log,
    Mul,
    Pow,
    preorder_traversal,
    sign,
    sin,
    simplify,
    srepr,
    S,
    sympify,
    var,
)
from nrpy.helpers.expr_tree import ExprTree
from nrpy.helpers.cse_preprocess_postprocess import cse_preprocess


# Basic Arithmetic Operations (Debugging)
def AbsSIMD_check(a: Symbol) -> Any:
    """Check that AbsSIMD(a) evaluates to sympy's Abs(a)."""
    return Abs(a)


def ConstSIMD_check(a: Basic) -> Float:
    """Convert a sympy Basic type to a Float with 34 digits of precision."""
    return Float(a, 34)


def nrpyAbsSIMD_check(a: Symbol) -> Any:
    """Check that nrpyAbsSIMD(a) evaluates to sympy's Abs(a)."""
    return Abs(a)


def AddSIMD_check(a: Symbol, b: Symbol) -> Any:
    """Perform addition operation for sympy Symbols."""
    return a + b


def SubSIMD_check(a: Symbol, b: Symbol) -> Any:
    """Perform subtraction operation for sympy Symbols."""
    return a - b


def MulSIMD_check(a: Symbol, b: Symbol) -> Any:
    """Perform multiplication operation for sympy Symbols."""
    return a * b


def FusedMulAddSIMD_check(a: Symbol, b: Symbol, c: Symbol) -> Any:
    """Perform fused multiply-add operation for sympy Symbols."""
    return a * b + c


def FusedMulSubSIMD_check(a: Symbol, b: Symbol, c: Symbol) -> Any:
    """Perform fused multiply-subtract operation for sympy Symbols."""
    return a * b - c


def NegFusedMulAddSIMD_check(a: Symbol, b: Symbol, c: Symbol) -> Any:
    """Perform negated fused multiply-add operation for sympy Symbols."""
    return -a * b + c


def NegFusedMulSubSIMD_check(a: Symbol, b: Symbol, c: Symbol) -> Any:
    """Perform negated fused multiply-subtract operation for sympy Symbols."""
    return -a * b - c


def DivSIMD_check(a: Symbol, b: Symbol) -> Any:
    """Perform division operation for sympy Symbols."""
    return a / b


def signSIMD_check(a: Basic) -> Any:
    """Find the sign of a sympy Basic type."""
    return sign(a)


# Transcendental Operations (Debugging)
def PowSIMD_check(a: Symbol, b: Symbol) -> Any:
    """Perform exponentiation operation for sympy Symbols."""
    return a**b


def SqrtSIMD_check(a: Symbol) -> Any:
    """Calculate the square root of a sympy Symbol."""
    return a ** (Rational(1, 2))


def CbrtSIMD_check(a: Symbol) -> Any:
    """Calculate the cube root of a sympy Symbol."""
    return a ** (Rational(1, 3))


def ExpSIMD_check(a: Symbol) -> Any:
    """Calculate the exponential of a sympy Symbol."""
    return exp(a)


def LogSIMD_check(a: Symbol) -> Any:
    """Calculate the natural logarithm of a sympy Symbol."""
    return log(a)


def SinSIMD_check(a: Symbol) -> Any:
    """Calculate the sine of a sympy Symbol."""
    return sin(a)


def CosSIMD_check(a: Symbol) -> Any:
    """Calculate the cosine of a sympy Symbol."""
    return cos(a)


def expr_convert_to_simd_intrins(
    expr: Basic,
    symbol_to_Rational_dict: Optional[Dict[Basic, Rational]] = None,
    prefix: str = "",
    simd_find_more_FMAsFMSs: bool = True,
    clean_NegativeOnes_after_processing: bool = False,
    debug: bool = False,
) -> Union[Basic, Expr]:
    """
    Convert expression to SIMD compiler intrinsics.

    :param expr: SymPy expression
    :param symbol_to_Rational_dict: Dictionary that maps symbols to rationals.
    :param prefix: Prefix for the function calls.
    :param simd_find_more_FMAsFMSs: Option to find more FMA/FMS patterns.
    :param clean_NegativeOnes_after_processing: Option to back-substitute and check difference.
    :param debug: Debug mode flag.
    :return: Expression containing SIMD compiler intrinsics.

    Doctests:
    >>> from sympy.abc import a, b, c, d
    >>> convert = expr_convert_to_simd_intrins

    >>> convert(-2*a)
    MulSIMD(-1, MulSIMD(_Integer_2, a))

    >>> convert(a**2)
    MulSIMD(a, a)

    >>> convert(a**(-2))
    DivSIMD(_Integer_1, MulSIMD(a, a))

    >>> convert(a**(1/2))
    SqrtSIMD(a)

    >>> convert(a**(-1/2))
    DivSIMD(_Integer_1, SqrtSIMD(a))

    >>> convert(a**(-3/2))
    DivSIMD(_Integer_1, MulSIMD(a, SqrtSIMD(a)))

    >>> convert(a**(-5/2))
    DivSIMD(_Integer_1, MulSIMD(MulSIMD(a, a), SqrtSIMD(a)))

    >>> convert(a**Rational(1, 3))
    CbrtSIMD(a)

    >>> convert(a**b)
    PowSIMD(a, b)

    >>> convert(a - b)
    SubSIMD(a, b)

    >>> convert(a + b - c)
    AddSIMD(b, SubSIMD(a, c))

    >>> convert(a + b + c)
    AddSIMD(a, AddSIMD(b, c))

    >>> convert(a + b + c + d)
    AddSIMD(AddSIMD(a, b), AddSIMD(c, d))

    >>> convert(a*b*c)
    MulSIMD(a, MulSIMD(b, c))

    >>> convert(a*b*c*d)
    MulSIMD(MulSIMD(a, b), MulSIMD(c, d))

    >>> convert(a/b)
    DivSIMD(a, b)

    >>> convert(a*b + c)
    FusedMulAddSIMD(a, b, c)

    >>> convert(a*b - c)
    FusedMulSubSIMD(a, b, c)

    >>> convert(-a*b + c)
    NegFusedMulAddSIMD(a, b, c)

    >>> convert(-a*b - c)
    NegFusedMulSubSIMD(a, b, c)
    """
    for item in preorder_traversal(expr):
        for arg in item.args:
            if isinstance(arg, Symbol):
                var(str(arg))

    if symbol_to_Rational_dict is None:
        cse_preprocessed_expr_list, symbol_to_Rational_dict = cse_preprocess(expr)
        expr = cse_preprocessed_expr_list[0]

    map_rat_to_sym = {symbol_to_Rational_dict[v]: v for v in symbol_to_Rational_dict}

    expr_orig, tree = expr, ExprTree(expr)

    AbsSIMD = Function("AbsSIMD")
    AddSIMD = Function("AddSIMD")
    SubSIMD = Function("SubSIMD")
    MulSIMD = Function("MulSIMD")
    FusedMulAddSIMD = Function("FusedMulAddSIMD")
    FusedMulSubSIMD = Function("FusedMulSubSIMD")
    NegFusedMulAddSIMD = Function("NegFusedMulAddSIMD")
    NegFusedMulSubSIMD = Function("NegFusedMulSubSIMD")
    DivSIMD = Function("DivSIMD")
    SignSIMD = Function("SignSIMD")

    PowSIMD = Function("PowSIMD")
    SqrtSIMD = Function("SqrtSIMD")
    CbrtSIMD = Function("CbrtSIMD")
    ExpSIMD = Function("ExpSIMD")
    LogSIMD = Function("LogSIMD")
    SinSIMD = Function("SinSIMD")
    CosSIMD = Function("CosSIMD")

    # Step 1: Replace transcendental functions, power functions, and division expressions.
    #   Note: SymPy does not represent fractional integers as rationals since
    #         those are explicitly declared using the rational class, and hence
    #         the following algorithm does not affect fractional integers.
    #         SymPy: srepr(a**(-2)) = Pow(a, -2)
    #         NRPy:  srepr(a**(-2)) = DivSIMD(1, MulSIMD(a, a))
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        if func == Abs:
            subtree.expr = AbsSIMD(args[0])
        elif func == exp:
            subtree.expr = ExpSIMD(args[0])
        elif func == log:
            subtree.expr = LogSIMD(args[0])
        elif func == sin:
            subtree.expr = SinSIMD(args[0])
        elif func == cos:
            subtree.expr = CosSIMD(args[0])
        elif func == sign:
            subtree.expr = SignSIMD(args[0])
    tree.reconstruct()

    def IntegerPowSIMD(a: Basic, n: int) -> Any:
        # Recursive Helper Function: Construct Integer Powers
        if n == 2:
            return MulSIMD(a, a)
        if n > 2:
            return MulSIMD(IntegerPowSIMD(a, n - 1), a)
        if n <= -2:
            one = Symbol(prefix + "_Integer_1")
            try:
                map_rat_to_sym[sympify(1)]
            except KeyError:
                symbol_to_Rational_dict[one], map_rat_to_sym[sympify(1)] = S.One, one
            return DivSIMD(one, IntegerPowSIMD(a, -n))
        if n == -1:
            one = Symbol(prefix + "_Integer_1")
            try:
                map_rat_to_sym[sympify(1)]
            except KeyError:
                symbol_to_Rational_dict[one], map_rat_to_sym[sympify(1)] = S.One, one
            return DivSIMD(one, a)
        if n == 1:
            return a
        # n == 0 case. Should never be called.
        return a

    def lookup_rational(
        arg: Union[Expr, Symbol, Rational, Basic]
    ) -> Union[Expr, Symbol, Rational, Basic]:
        if arg.func == Symbol:
            try:
                arg = symbol_to_Rational_dict[arg]
            except KeyError:
                pass
        return arg

    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        if func == Pow:
            one = Symbol(prefix + "_Integer_1")
            exponent = lookup_rational(args[1])
            if exponent == 0.5:
                subtree.expr = SqrtSIMD(args[0])
                subtree.children.pop(1)
            elif exponent == -0.5:
                subtree.expr = DivSIMD(one, SqrtSIMD(args[0]))
                tree.build(subtree)
            elif exponent == -1.5:
                pow_1p5 = SqrtSIMD(args[0]) * (args[0])
                subtree.expr = DivSIMD(one, pow_1p5)
                tree.build(subtree)
            elif exponent == -2.5:
                pow_2p5 = SqrtSIMD(args[0]) * (args[0]) * (args[0])
                subtree.expr = DivSIMD(one, pow_2p5)
                tree.build(subtree)
            elif exponent == Rational(1, 3):
                subtree.expr = CbrtSIMD(args[0])
                subtree.children.pop(1)
            elif isinstance(exponent, (Integer, int)):
                subtree.expr = IntegerPowSIMD(args[0], exponent)
                tree.build(subtree)
            else:
                subtree.expr = PowSIMD(*args)
    tree.reconstruct()

    # Step 2: Replace subtraction expressions.
    #   Note: SymPy: srepr(a - b) = Add(a, Mul(-1, b))
    #         NRPy:  srepr(a - b) = SubSIMD(a, b)
    for subtree in tree.preorder():
        func = subtree.expr.func
        args_list = list(subtree.expr.args)
        if func == Add:
            try:
                # Find the first occurrence of a negative product inside the addition
                i = next(
                    i
                    for i, arg in enumerate(args_list)
                    if arg.func == Mul
                    and any(lookup_rational(arg) == -1 for arg in args_list[i].args)
                )
                # Find the first occurrence of a negative symbol inside the product
                j = next(
                    j
                    for j, arg in enumerate(args_list[i].args)
                    if lookup_rational(arg) == -1
                )
                # Find the first non-negative argument of the product
                k = next(k for k in range(len(args_list)) if k != i)
                # Remove the negative symbol from the product
                subargs = list(args_list[i].args)
                subargs.pop(j)
                # Build the subtraction expression for replacement
                subexpr = SubSIMD(args_list[k], Mul(*subargs))
                args_list = [
                    arg for arg in args_list if arg not in (args_list[i], args_list[k])
                ]
                if len(args_list) > 0:
                    subexpr = Add(subexpr, *args_list)
                subtree.expr = subexpr
                tree.build(subtree)
            except StopIteration:
                pass
    tree.reconstruct()

    # Step 3: Replace addition and multiplication expressions.
    #   Note: SIMD addition and multiplication intrinsics can read
    #         only two arguments at once, whereas SymPy's Mul() and Add()
    #         operators can read an arbitrary number of arguments.
    #         SymPy: srepr(a*b*c*d) = Mul(a, b, d)
    #         NRPy:  srepr(a*b*c*d) = MulSIMD(MulSIMD(a, b), MulSIMD(c, d))
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        if func in (Mul, Add):
            func = MulSIMD if func == Mul else AddSIMD
            subexpr = func(*args[-2:])
            args, N = args[:-2], len(args) - 2
            for i in range(0, N, 2):
                if N - i > 1:
                    tmpexpr = func(args[i], args[i + 1])
                    subexpr = func(tmpexpr, subexpr, evaluate=False)
                else:
                    subexpr = func(args[i], subexpr, evaluate=False)
            subtree.expr = subexpr
            tree.build(subtree)
    tree.reconstruct()

    # Step 4: Replace the pattern Mul(Div(1, b), a) or Mul(a, Div(1, b)) with Div(a, b).
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # MulSIMD(DivSIMD(1, b), a) >> DivSIMD(a, b)
        if (
            func == MulSIMD
            and args[0].func == DivSIMD
            and lookup_rational(args[0].args[0]) == 1
        ):
            subtree.expr = DivSIMD(args[1], args[0].args[1])
            tree.build(subtree)
        # MulSIMD(a, DivSIMD(1, b)) >> DivSIMD(a, b)
        elif (
            func == MulSIMD
            and args[1].func == DivSIMD
            and lookup_rational(args[1].args[0]) == 1
        ):
            subtree.expr = DivSIMD(args[0], args[1].args[1])
            tree.build(subtree)
    tree.reconstruct()

    # Step 5: Now that all multiplication and addition functions only take two
    #         arguments, we can define fused-multiply-add functions,
    #         where AddSIMD(a, MulSIMD(b, c)) = b*c + a = FusedMulAddSIMD(b, c, a),
    #         or    AddSIMD(MulSIMD(b, c), a) = b*c + a = FusedMulAddSIMD(b, c, a).
    #   Note: Fused-multiply-add (FMA3) is standard on Intel CPUs with the AVX2
    #         instruction set, starting with Haswell processors in 2013:
    #         https://en.wikipedia.org/wiki/Haswell_(microarchitecture)

    # Step 5.a: Find double FMA patterns first [e.g. FMA(a, b, FMA(c, d, e))].
    #   Note: Double FMA simplifications do not guarantee a significant performance impact when solving BSSN equations
    if simd_find_more_FMAsFMSs:
        for subtree in tree.preorder():
            func = subtree.expr.func
            args = subtree.expr.args
            # a + b*c + d*e -> FMA(b,c,FMA(d,e,a))
            # AddSIMD(a, AddSIMD(MulSIMD(b,c), MulSIMD(d,e))) >> FusedMulAddSIMD(b, c, FusedMulAddSIMD(d,e,a))
            # Validate:
            # x = a + b*c + d*e
            # c_codegen(x,"x", params="enable_simd=True,SIMD_debug=True")
            if (
                func == AddSIMD
                and args[1].func == AddSIMD
                and args[1].args[0].func == MulSIMD
                and args[1].args[1].func == MulSIMD
            ):
                subtree.expr = FusedMulAddSIMD(
                    args[1].args[0].args[0],
                    args[1].args[0].args[1],
                    FusedMulAddSIMD(
                        args[1].args[1].args[0], args[1].args[1].args[1], args[0]
                    ),
                )
                tree.build(subtree)
            # b*c + d*e + a -> FMA(b,c,FMA(d,e,a))
            # Validate:
            # x = b*c + d*e + a
            # c_codegen(x,"x", params="enable_simd=True,SIMD_debug=True")
            # AddSIMD(AddSIMD(MulSIMD(b,c), MulSIMD(d,e)),a) >> FusedMulAddSIMD(b, c, FusedMulAddSIMD(d,e,a))
            elif (
                func == AddSIMD
                and args[0].func == AddSIMD
                and args[0].args[0].func == MulSIMD
                and args[0].args[1].func == MulSIMD
            ):
                subtree.expr = FusedMulAddSIMD(
                    args[0].args[0].args[0],
                    args[0].args[0].args[1],
                    FusedMulAddSIMD(
                        args[0].args[1].args[0], args[0].args[1].args[1], args[1]
                    ),
                )
                tree.build(subtree)
        tree.reconstruct()

    # Step 5.b: Find single FMA patterns.
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # AddSIMD(MulSIMD(b, c), a) >> FusedMulAddSIMD(b, c, a)
        if func == AddSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree)
        # AddSIMD(a, MulSIMD(b, c)) >> FusedMulAddSIMD(b, c, a)
        elif func == AddSIMD and args[1].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[1].args[0], args[1].args[1], args[0])
            tree.build(subtree)
        # SubSIMD(MulSIMD(b, c), a) >> FusedMulSubSIMD(b, c, a)
        elif func == SubSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulSubSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree)
        # SubSIMD(a, MulSIMD(b, c)) >> NegativeFusedMulAddSIMD(b, c, a)
        elif func == SubSIMD and args[1].func == MulSIMD:
            subtree.expr = NegFusedMulAddSIMD(args[1].args[0], args[1].args[1], args[0])
            tree.build(subtree)
        # FMS(-1, MulSIMD(a, b), c) >> NegativeFusedMulSubSIMD(b, c, a)
        func = subtree.expr.func
        args = subtree.expr.args
        if (
            func == FusedMulSubSIMD
            and args[1].func == MulSIMD
            and lookup_rational(args[0]) == -1
        ):
            subtree.expr = NegFusedMulSubSIMD(args[1].args[0], args[1].args[1], args[2])
            tree.build(subtree)
    tree.reconstruct()

    # Step 5.c: Remaining double FMA patterns that previously in Step 5.a were difficult to find.
    #   Note: Double FMA simplifications do not guarantee a significant performance impact when solving BSSN equations
    if simd_find_more_FMAsFMSs:
        for subtree in tree.preorder():
            func = subtree.expr.func
            args = subtree.expr.args
            # (b*c - d*e) + a -> AddSIMD(a, FusedMulSubSIMD(b, c, MulSIMD(d, e))) >> FusedMulSubSIMD(b, c, FusedMulSubSIMD(d,e,a))
            # Validate:
            # x = (b*c - d*e) + a
            # c_codegen(x,"x", params="enable_simd=True,SIMD_debug=True")
            if (
                func == AddSIMD
                and args[1].func == FusedMulSubSIMD
                and args[1].args[2].func == MulSIMD
            ):
                subtree.expr = FusedMulSubSIMD(
                    args[1].args[0],
                    args[1].args[1],
                    FusedMulSubSIMD(
                        args[1].args[2].args[0], args[1].args[2].args[1], args[0]
                    ),
                )
                tree.build(subtree)
            # b*c - (a - d*e) -> SubSIMD(FusedMulAddSIMD(b, c, MulSIMD(d, e)), a) >> FMA(b,c,FMS(d,e,a))
            # Validate:
            # x = b * c - (a - d * e)
            # c_codegen(x, "x", params="enable_simd=True,SIMD_debug=True")
            elif (
                func == SubSIMD
                and args[0].func == FusedMulAddSIMD
                and args[0].args[2].func == MulSIMD
            ):
                subtree.expr = FusedMulAddSIMD(
                    args[0].args[0],
                    args[0].args[1],
                    FusedMulSubSIMD(
                        args[0].args[2].args[0], args[0].args[2].args[1], args[1]
                    ),
                )
                tree.build(subtree)
            # (b*c - d*e) - a -> SubSIMD(FusedMulSubSIMD(b, c, MulSIMD(d, e)), a) >> FMS(b,c,FMA(d,e,a))
            # Validate:
            # x = (b*c - d*e) - a
            # c_codegen(x,"x", params="enable_simd=True,SIMD_debug=True")
            elif (
                func == SubSIMD
                and args[0].func == FusedMulSubSIMD
                and args[0].args[2].func == MulSIMD
            ):
                subtree.expr = FusedMulSubSIMD(
                    args[0].args[0],
                    args[0].args[1],
                    FusedMulAddSIMD(
                        args[0].args[2].args[0], args[0].args[2].args[1], args[1]
                    ),
                )
                tree.build(subtree)
        tree.reconstruct()

    # Step 5.d: NegFusedMulAddSIMD(a,b,c) = -a*b + c:
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # FMA(a,Mul(-1,b),c) >> NFMA(a,b,c)
        if (
            func == FusedMulAddSIMD
            and args[1].func == MulSIMD
            and lookup_rational(args[1].args[0]) == -1
        ):
            subtree.expr = NegFusedMulAddSIMD(args[0], args[1].args[1], args[2])
            tree.build(subtree)
        # FMA(a,Mul(b,-1),c) >> NFMA(a,b,c)
        elif (
            func == FusedMulAddSIMD
            and args[1].func == MulSIMD
            and lookup_rational(args[1].args[1]) == -1
        ):
            subtree.expr = NegFusedMulAddSIMD(args[0], args[1].args[0], args[2])
            tree.build(subtree)
        # FMA(Mul(-1,a), b,c) >> NFMA(a,b,c)
        elif (
            func == FusedMulAddSIMD
            and args[0].func == MulSIMD
            and lookup_rational(args[0].args[0]) == -1
        ):
            subtree.expr = NegFusedMulAddSIMD(args[0].args[1], args[1], args[2])
            tree.build(subtree)
        # FMA(Mul(a,-1), b,c) >> NFMA(a,b,c)
        elif (
            func == FusedMulAddSIMD
            and args[0].func == MulSIMD
            and lookup_rational(args[0].args[1]) == -1
        ):
            subtree.expr = NegFusedMulAddSIMD(args[0].args[0], args[1], args[2])
            tree.build(subtree)
    tree.reconstruct()

    # Step 5.e: Replace e.g., FMA(-1,b,c) with SubSIMD(c,b) and similar patterns
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # FMA(-1,b,c) >> SubSIMD(c,b)
        if func == FusedMulAddSIMD and lookup_rational(args[0]) == -1:
            subtree.expr = SubSIMD(args[2], args[1])
            tree.build(subtree)
        # FMA(a,-1,c) >> SubSIMD(c,a)
        elif func == FusedMulAddSIMD and lookup_rational(args[1]) == -1:
            subtree.expr = SubSIMD(args[2], args[0])
            tree.build(subtree)
        # FMS(a,-1,c) >> MulSIMD(-1,AddSIMD(a,c))
        elif func == FusedMulSubSIMD and lookup_rational(args[1]) == -1:
            subtree.expr = MulSIMD(args[1], AddSIMD(args[0], args[2]))
            tree.build(subtree)
        # FMS(-1,b,c) >> MulSIMD(-1,AddSIMD(b,c))
        elif func == FusedMulSubSIMD and lookup_rational(args[0]) == -1:
            subtree.expr = MulSIMD(args[0], AddSIMD(args[1], args[2]))
            tree.build(subtree)
    tree.reconstruct()

    # Step 5.f: NegFusedMulSubSIMD(a,b,c) = -a*b - c:
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # NFMA(a,b,Mul(-1,c)) >> NFMS(a,b,c)
        if (
            func == NegFusedMulAddSIMD
            and args[2].func == MulSIMD
            and lookup_rational(args[2].args[0]) == -1
        ):
            subtree.expr = NegFusedMulSubSIMD(args[0], args[1], args[2].args[1])
            tree.build(subtree)
        # NFMA(a,b,Mul(c,-1)) >> NFMS(a,b,c)
        elif (
            func == NegFusedMulAddSIMD
            and args[2].func == MulSIMD
            and lookup_rational(args[2].args[1]) == -1
        ):
            subtree.expr = NegFusedMulSubSIMD(args[0], args[1], args[2].args[0])
            tree.build(subtree)
        # FMS(a,Mul(-1,b),c) >> NFMS(a,b,c)
        elif (
            func == FusedMulSubSIMD
            and args[1].func == MulSIMD
            and lookup_rational(args[1].args[0]) == -1
        ):
            subtree.expr = NegFusedMulSubSIMD(args[0], args[1].args[1], args[2])
            tree.build(subtree)
        # FMS(a,Mul(b,-1),c) >> NFMS(a,b,c)
        elif (
            func == FusedMulSubSIMD
            and args[1].func == MulSIMD
            and lookup_rational(args[1].args[1]) == -1
        ):
            subtree.expr = NegFusedMulSubSIMD(args[0], args[1].args[0], args[2])
            tree.build(subtree)
        # FMS(a,Mul([something],Mul(-1,b)),c) >> NFMS(a,Mul([something],b),c)
        elif (
            func == FusedMulSubSIMD
            and args[1].func == MulSIMD
            and args[1].args[1].func == MulSIMD
            and lookup_rational(args[1].args[1].args[0]) == -1
        ):
            subtree.expr = NegFusedMulSubSIMD(
                args[0], MulSIMD(args[1].args[0], args[1].args[1].args[1]), args[2]
            )
            tree.build(subtree)
        # FMS(a,Mul([something],Mul(b,-1)),c) >> NFMS(a,Mul([something],b),c)
        elif (
            func == FusedMulSubSIMD
            and args[1].func == MulSIMD
            and args[1].args[1].func == MulSIMD
            and lookup_rational(args[1].args[1].args[1]) == -1
        ):
            subtree.expr = NegFusedMulSubSIMD(
                args[0], MulSIMD(args[1].args[0], args[1].args[1].args[0]), args[2]
            )
            tree.build(subtree)
    tree.reconstruct()

    # Step 5.g: Find single FMA patterns again, as some new ones might be found.
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # AddSIMD(MulSIMD(b, c), a) >> FusedMulAddSIMD(b, c, a)
        if func == AddSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree)
        # AddSIMD(a, MulSIMD(b, c)) >> FusedMulAddSIMD(b, c, a)
        elif func == AddSIMD and args[1].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[1].args[0], args[1].args[1], args[0])
            tree.build(subtree)
        # SubSIMD(MulSIMD(b, c), a) >> FusedMulSubSIMD(b, c, a)
        elif func == SubSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulSubSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree)
    expr = tree.reconstruct()

    if debug:
        # pylint: disable=W0123
        expr_check = eval(str(expr).replace("SIMD", "SIMD_check"))
        expr_check = expr_check.subs(-1, Symbol("_NegativeOne_"))

        expr_diff = expr_check - expr_orig
        # The eval(str(srepr())) below normalizes the expression,
        # fixing a cancellation issue in SymPy ~0.7.4.
        # pylint: disable=W0123
        expr_diff = eval(str(srepr(expr_diff)))
        tree_diff = ExprTree(expr_diff)
        for subtree in tree_diff.preorder():
            subexpr = subtree.expr
            if subexpr.func == Float:
                if abs(subexpr - Integer(subexpr)) < 1.0e-14 * subexpr:
                    subtree.expr = Integer(subexpr)
        expr_diff = tree_diff.reconstruct()

        if expr_diff != 0:
            # MyPy gives a completely bogus error: Module not callable  [operator], referring to simplify().
            simp_expr_diff = simplify(expr_diff)
            if simp_expr_diff != 0:
                raise Warning("Expression Difference: " + str(simp_expr_diff))

    def remove_unused_NegativeOnes_from_symbol_to_Rational_dict(
        expr: Union[Basic, Expr], symbol_to_Rational_dict: Dict[Basic, Rational]
    ) -> None:
        """
        In matching many patterns above, we have removed NegativeOne's from expressions.
        If all -1 have been removed, this function removes {prefix}_NegativeOne_ from
        symbol_to_Rational_dict, so it doesn't get declared as an unused temporary variable.

        :param expr: The mathematical expression from which NegativeOne's have been removed.
        :param symbol_to_Rational_dict: Dictionary mapping symbols to Rational numbers.
        """
        NegOne_in_symbol_to_Rational_dict = False
        NegOne_symb = Symbol("none") * 2
        NegOne_symb_str = f"{prefix}_NegativeOne_"
        for symb in symbol_to_Rational_dict.keys():
            if str(symb) == NegOne_symb_str:
                NegOne_symb = symb
                NegOne_in_symbol_to_Rational_dict = True
                break
        found_NegOne_in_free_symbols = False
        for symb in expr.free_symbols:
            if str(symb) == NegOne_symb_str:
                found_NegOne_in_free_symbols = True
                break
        if NegOne_in_symbol_to_Rational_dict and not found_NegOne_in_free_symbols:
            del symbol_to_Rational_dict[NegOne_symb]

    if symbol_to_Rational_dict and clean_NegativeOnes_after_processing:
        remove_unused_NegativeOnes_from_symbol_to_Rational_dict(
            expr, symbol_to_Rational_dict
        )

    return expr


def copy_simd_intrinsics_h(project_dir: str) -> None:
    """
    Copy simd_intrinsics.h into the specified project directory.

    :param project_dir: The path of the project directory where the file will be copied.
    :return: None

    :raises ImportError: If necessary libraries are not found.
    :raises FileNotFoundError: If source file is not found.
    """
    simd_path = Path(project_dir) / "simd"
    simd_path.mkdir(parents=True, exist_ok=True)

    try:
        # only Python 3.7+ has importlib.resources
        import importlib.resources  # pylint: disable=E1101,C0415

        source_path = (
            # pylint: disable=E1101
            importlib.resources.files("nrpy.helpers")
            / "simd_intrinsics.h"
        )
        shutil.copy(str(source_path), str(simd_path))
    except ImportError:  # Fallback to resource_filename for older Python versions
        # pylint: disable=E1101,C0415
        from pkg_resources import resource_filename  # type: ignore

        source_path = resource_filename("nrpy.helpers", "simd_intrinsics.h")
        shutil.copy(source_path, str(simd_path))  # type: ignore


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
