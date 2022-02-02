# py-double-exponential
A quick and dirty double-exponential (aka `tanh-sinh`) integration program in Python using `mpmath`.

There are 3 Python files here:

- `double_exponential.py`: this contains the (quick) function that performs the quadrature. Read its docstring for usage tips. This file can also be used from the command line (see below)
- `double_exponential_tests.py`: this (dirty) script uses the previous function to evaluate some test integrals and report the achieved results (see the docstring for a description of the output format). The test integrals are defined in…
- `test_integrals.py`: contains a list of use cases to test the algorithm. Read its docstring to get the format in order to add more use cases.

Uses `mpmath` and the only adjustable parameter is the number of bits or decimal digits used during calculations. It (`mp.mp.prec` or `mp.mp.dps`) can be adjusted at the beginning of `double_exponential_tests.py` or before entering `double_exponential.py`.

This algorithm accepts `±mp.inf` as integration limits.

## Using `double_exponential.py` from the command line
`double_exponential.py` can be invoked from the command line. Its usage is:

```
  double_exponential.py [-h] [-b BITS] [-d DIGITS] f a b
```

with positional arguments:
  * `f`: function to be integrated, may be a " quoted string with a lambda expression, e.g.: "`lambda x: 2/(1 + x**2)`"
  * `a`: left extreme of the integration interval (accepts infinities as "`±mp.inf`")
  * `b`: right extreme of the integration interval (accepts infinities as "`±mp.inf`")

and options:
  * `-h`, `--help`: shows a help message and exits
  * `-b BITS`, `--bits BITS`: sets the number of bits used during calculations (sets `mp.prec`)
  * `-d DIGITS`, `--digits DIGITS`: sets the number of decimal digits used during calculations (sets `mp.dps`)

At end it will return the computed integral, the estimated error of the quadrature, the Total Number of Function Evaluations (TNFE)
needed during the calculation and the used variant of the method. When specifying the arguments, it can be assumed
that package `mpmath` is imported as `mp`. Use "`--`" (once) after that last option and before any argument that starts with a "`-`". Thus:
```
  double_exponential.py -d 32 -- "lambda x: 2/(1 + x**2)" -mp.inf 0
```
will return:
```
I = 3.1415926535897932384626433832795 ± 9.1769596347412488797710510388239e-16

TNFE = 95 (exp-sinh)
```
Note that, in this case, all 32 digits of the result are correct.
