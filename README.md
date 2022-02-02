# py-double-exponential

A quick and dirty double-exponential (aka `tanh-sinh`) integration program in Python.

There are 3 Python files here:

- `double_exponential.py`: this contains the (quick) function that performs the quadrature. Read its docstring for usage tips.
- `double_exponential_tests.py`: this (dirty) script uses the previous function to evaluate some test integrals and report the achieved results (see the docstring for a description of the output format). The test integrals are defined inâ€¦
- `test_integrals.py`: contains a list of use cases to test the algorithm. Read its docstring to get the format in order to add more use cases.

Uses `mpmath` and the only adjustable parameter is the number of decimal digits used during calculations. It (`mp.dps`) can be adjusted at the beginning of `double_exponential_tests.py`.

This algorithm accepts `+inf`/`-inf` as integration limits.

Enjoy!

