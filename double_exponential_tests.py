# Copyright (c) 2021, emece67 - MIT License


"""
Transverses a list `test_integral` of test-cases computing the specified
  integral calling double_exponential(). For each case reports:
    * #:      case number
    * TNFE:   Total Number of Function Evaluations
    * Lvl:    level reached during iteration
    * CD:     correct digits achieved
    * err_r:  reported error
    * err_t:  true error
    * q:      computed integral
    * I:      true integral

At end reports accumulated results:
    * TNFE
    * CD
    * ratio of CD to total digits
  for the 3 variants of the algorithm (sinh-sinh, exp-sinh and tanh-sinh).

Set mp.dps to the desired number of decimal digits to be used in
  calculations.
"""


from mpmath import mp, isfinite, floor, log10


# number of decimal digits used in computations
mp.dps = 15


from double_exponential import double_exponential
from test_integrals import test_integral


# #:      Number
# TNFE:   Total Number of Function Evaluations
# Lvl:    level reached during iteration
# CD:     Correct Digits
# err_r:  reported error
# err_t:  true error
# q:      computed integral
# I:      true integral
print('#   TNFE Lvl CD err_r err_t %s I'%('q'.ljust(mp.dps + 6)))

# accumulated results
res = {
  'SS_TNFE' : 0,    # accumulated TNFE for all sinh-sinh cases
  'SS_CD': 0,       # accumulated CD for all sinh-sinh cases
  'SS_n': 0,        # number of sinh-sinh cases
  'ES_TNFE': 0,     # same for exp-sinh
  'ES_CD': 0,
  'ES_n': 0,
  'TS_TNFE': 0,     # and for tanh-sinh
  'TS_CD': 0,
  'TS_n': 0}

# try (some/all) test cases
for n, integral in enumerate(test_integral[:]):
  q, err_r, tnfe, lvl, variant, q_lvl = double_exponential(integral['f'], \
                                          integral['a'], integral['b'])
  I = integral['s']

  # number of correct digits (CD)
  if I != 0 and isfinite(I):
    cd =  max(\
              min(int(floor(log10(abs(I))) - floor(log10(abs(q - I)))), \
                  mp.dps), \
              0) if q != I \
          else mp.dps
  elif I == 0:
    cd = int(-round(log10(abs(q)))) if q != 0 else mp.dps
  else:
    cd = 0

  # populate results
  if variant == 0:
    res['TS_TNFE'] += tnfe
    res['TS_CD'] += cd
    res['TS_n'] += 1
  elif variant == 1:
    res['ES_TNFE'] += tnfe
    res['ES_CD'] += cd
    res['ES_n'] += 1
  else:
    res['SS_TNFE'] += tnfe
    res['SS_CD'] += cd
    res['SS_n'] += 1

  # print results for each integral
  if I != 0 and isfinite(I):
    print('%03i %04i %02i '%(n, tnfe, lvl), end = ' ')
    print('%02i'%(cd), end = ' ')
    print('%1.0e'%(err_r), end = ' ')
    print('%1.0e'%(abs(I - q)), end = ' ')
    print('%+*.*e'%(mp.dps + 6, mp.dps - 1, q), \
            '%+*.*e'%(mp.dps + 6, mp.dps - 1, I), end = '')
  elif I == 0:
    print('%03i %04i %02i  %02i %1.0e %01.0e'%(n, tnfe, lvl, cd, err_r, abs(I - q)), \
            '%+*.*e'%(mp.dps + 6, mp.dps - 1, q), '0', end = '')
  else:
    print('%03i %04i %02i  %02i inf   inf  '%(n, tnfe, lvl, cd), \
            '%+*.*e'%(mp.dps + 6, mp.dps - 1, q), 'inf', end = '')
  print()

# print cumulative results
res['SS_n'] += 1 if not res['SS_n'] else 0
res['ES_n'] += 1 if not res['ES_n'] else 0
res['TS_n'] += 1 if not res['TS_n'] else 0
print()
print('SS: %05i - %04i (%02i)\nES: %05i - %04i (%02i)\nTS: %05i - %04i (%02i)\n'
  %(res['SS_TNFE'], res['SS_CD'], 100*res['SS_CD']/mp.dps/res['SS_n'],
    res['ES_TNFE'], res['ES_CD'], 100*res['ES_CD']/mp.dps/res['ES_n'],
    res['TS_TNFE'], res['TS_CD'], 100*res['TS_CD']/mp.dps/res['TS_n']))
