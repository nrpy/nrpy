"""
Functional Programming Toolkit.

Author: Ken Sible
Email:  ksible *at* outlook *dot* com
"""

import sys
from typing import Callable, Iterable, List, Tuple, Optional, Any, Generator, Union


def pipe(x: Any, *f: Callable[..., Any]) -> Any:
    """
    Pipe operator that applies a sequence of functions to an input.

    :param x: The input to the first function in the sequence.
    :param f: Functions to apply in sequence.
    :return: The result of applying all the functions to the input.

    Doctest:
    >>> pipe(range(5, 0, -1), reversed, list)
    [1, 2, 3, 4, 5]
    >>> pipe([3, 2, 2, 4, 5, 1], sorted, set, list)
    [1, 2, 3, 4, 5]
    """
    if not f:
        return x
    return pipe(f[0](x), *f[1:])


def repeat(f: Callable[..., Any], x: Any, n: int) -> Any:
    """
    Repeats the application of a function a given number of times.

    :param f: The function to apply.
    :param x: The input to the function.
    :param n: The number of times to apply the function.
    :return: The result of repeatedly applying the function.

    Doctest:
    >>> list(repeat(flatten, [1, 2, [3, [4]], 5], 2))
    [1, 2, 3, 4, 5]
    """
    if n == 0:
        return x
    return repeat(f, f(x), n - 1)


def chain(*iterable: Iterable[Any]) -> Generator[Any, None, None]:
    """
    Chains several iterables together.

    :param iterable: A sequence of iterables to chain together.
    :return: A generator that yields elements from each iterable in sequence.

    Doctest:
    >>> list(chain([1], [2, 3], [4, 5]))
    [1, 2, 3, 4, 5]
    """
    for iter_ in iterable:
        try:
            iter(iter_)
        except TypeError:
            iter_ = [iter_]
        for element in iter_:
            yield element


def flatten(iterable: Iterable[Any]) -> Generator[Any, None, None]:
    """
    Flattens a nested iterable into a single-level iterable.

    :param iterable: The iterable to flatten.
    :return: A generator that yields the flattened elements of the iterable.

    Doctest:
    >>> list(flatten([1, [2, 3], [4, 5]]))
    [1, 2, 3, 4, 5]
    """
    return chain(*iterable)


def reduce(
    f: Callable[[Any, Any], Any],
    iterable: Iterable[Any],
    initializer: Optional[Any] = None,
) -> Any:
    """
    Apply a binary function cumulatively to all the elements of an iterable.

    :param f: The binary function to apply.
    :param iterable: The iterable to reduce.
    :param initializer: The initial value to start the reduction.
    :return: The result of reducing the iterable.

    Doctests:
    >>> reduce(lambda x, y: x + y, [1, 2, 3, 4, 5])
    15
    >>> reduce(lambda x, y: x + y, ['w', 'o', 'r', 'd'])
    'word'
    >>> x = [1, 2, [3, 4], 'aabb']
    >>> reduce(lambda i, _: i + 1, [1] + x[1:])
    4
    """
    iterable = iter(iterable)
    result = next(iterable) if initializer is None else initializer
    for element in iterable:
        result = f(result, element)
    return result


def uniquify(iterable: Iterable[Any]) -> List[Any]:
    """
    Return a list with the unique elements of the input iterable, maintaining the order.

    :param iterable: An iterable.
    :return: A list with the unique elements of the iterable.

    Doctest:
    >>> uniquify(([1, 1, 2, 3, 3, 3, 4, 5, 5]))
    [1, 2, 3, 4, 5]
    """
    return list(reduce(lambda l, x: l if x in l else l + [x], iterable, []))


def product(
    *iterable: Union[Iterable[Any], Tuple[Iterable[Any], ...]], **kwargs: Any
) -> Generator[Tuple[Any, ...], None, None]:
    """
    Return the Cartesian product of several iterables.

    :param iterable: A sequence of iterables.
    :param kwargs: An optional keyword argument 'repeat' to specify a repetition of a single iterable.
    :return: A generator that yields tuples, each containing one element from each input iterable.

    Doctests:
    >>> list(product(['a', 'b'], [1, 2, 3]))
    [('a', 1), ('a', 2), ('a', 3), ('b', 1), ('b', 2), ('b', 3)]
    >>> list(product([1, 2, 3], repeat=2))
    [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
    >>> for i, j in product(['a', 'b'], range(1, 3)):
    ...     print('%s: %d' % (i, j))
    ...
    a: 1
    a: 2
    b: 1
    b: 2
    """
    if "repeat" in kwargs:
        if kwargs["repeat"] > 1 and len(iterable) == 1:
            iterable = kwargs["repeat"] * iterable

    def f(A: List[List[Any]], B: List[List[Any]]) -> List[List[Any]]:
        return [list(flatten([x] + [y])) for x in A for y in B]

    for prod in reduce(f, iterable):
        yield tuple(prod)


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
