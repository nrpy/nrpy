"""
Symbolic (N-Ary) Expression Tree.

The following script extends the expression tree from SymPy,
allowing direct node manipulation for subexpression replacement.
The expression tree structure within SymPy expressions stores
subexpressions inside immutable tuples, which prevents clients from
modifying the expression tree. Therefore, clients must depend on
built-in functions, such as xreplace, for subexpression replacement,
which might be suboptimal for their specific purposes. The ExprTree class
is implemented as an n-ary tree data structure for SymPy expressions,
and is equipped with a build method for constructing the expression tree,
a reconstruct method for reconstructing the root expression, a replace
method for subexpression replacement, and preorder/postorder traversal
iterators (or generators). The __repr__ method of the expression
tree returns a string of the expressions using the preorder traversal,
while the __str__ method returns a string of the class name
and root expression. The Node subclass has a field for an expression and
a field for subexpression children, which is implemented as a mutable list.

Author: Ken Sible
Email:  ksible *at* outlook *dot* com
"""

import sys  # Standard Python module for multiplatform OS-level functions
from typing import Generator, List, Optional, Union

import sympy as sp


class ExprTree:
    """
    Represents a symbolic expression tree for a given SymPy expression.

    This class is designed to facilitate tree-based manipulation and traversal
    of symbolic expressions, supporting operations like building the tree from
    an expression, preorder and postorder traversal.

    Attributes
    ----------
    root : Node
        The root node of the expression tree.

    Example:
    >>> from sympy.abc import a, b
    >>> from sympy import cos
    >>> tree = ExprTree(cos(a + b)**2)
    >>> print(tree)
    ExprTree(cos(a + b)**2)
    >>> [node.expr for node in tree.preorder()]
    [cos(a + b)**2, cos(a + b), a + b, a, b, 2]
    """

    class Node:
        """
        Represents a node within an `ExprTree`, holding a SymPy expression and its function.

        Attributes
        ----------
        expr : sp.Basic
            The SymPy expression held by this node.
        func : Optional[sp.FunctionClass]
            The function of the expression if applicable.
        children : List["ExprTree.Node"]
            Child nodes under this node.
        """

        def __init__(self, expr: sp.Basic, func: Optional[sp.FunctionClass]) -> None:
            """
            Initialize an `ExprTree.Node` instance.

            :param expr: The SymPy expression for this node.
            :param func: The function of the expression, if applicable.
            """
            self.expr: sp.Basic = expr
            self.func = func
            self.children: List["ExprTree.Node"] = []

        def append(self, node: "ExprTree.Node") -> None:
            """
            Append a child node to this node.

            :param node: The `ExprTree.Node` instance to append as a child.
            """
            self.children.append(node)

    def __init__(self, expr: sp.Basic) -> None:
        """
        Initialize an `ExprTree` with a given root expression.

        :param expr: The SymPy expression to set as the root of the tree.
        """
        self.root = self.Node(expr, None)
        self.build(self.root)

    def build(self, node: "ExprTree.Node", clear: bool = True) -> None:
        """
        Recursively build the expression tree from a given node.

        :param node: The node to start building the tree from.
        :param clear: If `True`, clears existing children of `node` before building.

        Example:
        >>> from sympy.abc import a, b
        >>> from sympy import cos, sin
        >>> tree = ExprTree(cos(a + b)**2)
        >>> tree.root.expr = sin(a*b)**2
        >>> tree.build(tree.root, clear=True)
        >>> [node.expr for node in tree.preorder()]
        [sin(a*b)**2, sin(a*b), a*b, a, b, 2]
        """
        if clear:
            del node.children[:]
        for arg in node.expr.args:
            subtree = self.Node(arg, node.expr.func)
            node.append(subtree)
            self.build(subtree)

    def preorder(
        self, node: Optional["ExprTree.Node"] = None
    ) -> Generator["ExprTree.Node", None, None]:
        """
        Perform a preorder traversal of the tree starting from a given node.

        :param node: The starting node for the traversal. If `None`, starts from the root.
        :yields: Each node in preorder sequence.

        Example:
        >>> from sympy.abc import a, b
        >>> from sympy import cos, Mul
        >>> tree = ExprTree(cos(a*b)**2)
        >>> for i, subtree in enumerate(tree.preorder()):
        ...     if subtree.expr.func == Mul:
        ...         print((i, subtree.expr))
        (2, a*b)
        """
        if node is None:
            node = self.root
        yield node
        for child in node.children:
            yield from self.preorder(child)

    def postorder(
        self, node: Optional["ExprTree.Node"] = None
    ) -> Generator["ExprTree.Node", None, None]:
        """
        Perform a postorder traversal of the tree starting from a given node.

        :param node: The starting node for the traversal. If `None`, starts from the root.
        :yields: Each node in postorder sequence.

        Example:
        >>> from sympy.abc import a, b
        >>> from sympy import cos, Mul
        >>> tree = ExprTree(cos(a*b)**2)
        >>> for i, subtree in enumerate(tree.postorder()):
        ...     if subtree.expr.func == Mul:
        ...         print((i, subtree.expr))
        (2, a*b)
        """
        if node is None:
            node = self.root
        for child in node.children:
            yield from self.postorder(child)
        yield node

    def reconstruct(self, evaluate: bool = False) -> sp.Basic:
        """
        Reconstruct root expression from expression tree.

        :param evaluate: Evaluate root expression, default is False
        :return: Root expression

        >>> from sympy.abc import a, b
        >>> from sympy import cos, sin
        >>> tree = ExprTree(cos(a + b)**2)
        >>> tree.root.children[0].expr = sin(a + b)
        >>> tree.reconstruct()
        sin(a + b)**2
        """
        for subtree in self.postorder():
            if subtree.children:
                expr_list = [node.expr for node in subtree.children]
                subtree.expr = subtree.expr.func(*expr_list, evaluate=evaluate)
        return self.root.expr

    def __repr__(self) -> str:
        return "ExprTree(" + str(self.root.expr) + ")"

    __str__ = __repr__

if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
