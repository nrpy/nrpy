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
from typing import Generator, List, Optional
import sympy as sp


class ExprTree:
    """
    Symbolic (N-Ary) Expression Tree.

    :param expr: Root expression for the tree

    >>> from sympy.abc import a, b, x
    >>> from sympy import cos
    >>> tree = ExprTree(cos(a + b)**2)
    >>> print(tree)
    ExprTree(cos(a + b)**2)
    >>> [node.expr for node in tree.preorder()]
    [cos(a + b)**2, cos(a + b), a + b, a, b, 2]
    """

    class Node:
        """
        Expression Tree Node.

        :param expr: Expression for the node
        :param func: Function for the expression
        """

        def __init__(
            self,
            expr: sp.Basic,
            func: Optional[sp.FunctionClass],
        ) -> None:
            self.expr: sp.Basic = expr
            self.func = func
            self.children: List["ExprTree.Node"] = []

        def append(self, node: "ExprTree.Node") -> None:
            """Append to a node in expression tree."""
            self.children.append(node)

        def __repr__(self) -> str:
            return f"Node({self.expr}, {self.func})"

        def __str__(self) -> str:
            return str(self.expr)

    def __init__(self, expr: sp.Basic) -> None:
        self.root = self.Node(expr, None)
        self.build(self.root)

    def build(self, node: Node, clear: bool = True) -> None:
        """
        Build expression (sub)tree.

        :param node: Root node of (sub)tree
        :param clear: Clear children, default is True

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

    def preorder(self, node: Optional[Node] = None) -> Generator[Node, None, None]:
        """
        Generate iterator for preorder traversal.

        :param node: Root node of (sub)tree
        :return: Iterator

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
            for subtree in self.preorder(child):
                yield subtree

    def postorder(self, node: Optional[Node] = None) -> Generator[Node, None, None]:
        """
        Generate iterator for postorder traversal.

        :param node: Root node of (sub)tree
        :return: Iterator

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
            for subtree in self.postorder(child):
                yield subtree
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
