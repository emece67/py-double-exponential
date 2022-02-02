# Copyright (c) 2021, emece67 - MIT License


from mpmath import isnan, isfinite, sign, almosteq, power, mp, sqrt, log, ln, pi, exp, isnormal


def double_exponential(f, a, b):
  """
  Computes the integral of function `f` from `a` to `b`, using the double
    exponential method. Accepts `+inf`/`-inf` as interval ends
  Returns a tuple with:
    * the computed integral;
    * an error estimation;
    * the number of function evaluations needed for the computation;
    * the level reached during computation;
    * the variant used: tanh-sinh (0); exp-sinh (1); or sinh-sinh (2)
    * a list with the approximations reached at each level

  Uses mpmath and works with the pre-existing `mp.dps` precision.

  If the computed error estimation is not much smaller than the computed
    result, it is assumed that all digits of the result are corrupted by
    roundoff. In such cases, the reported result is 0 and the reported
    error estimation equals the sum of the absolute values of the
    computed result and error. This usually happens when the integral
    evaluates to 0.

  The double exponential method relies on the function to be
    integrated being analytic over the integration interval, except,
    perhaps, at the interval ends. Thus, if it is known that the
    integrand, _or any of its derivatives_, has a discontinuity at some
    point inside the integration interval, it is advised to break the
    interval at such point and compute two separate integrals, one
    at each side of the problematic point. Better results will be
    get this way. Thus, beware of `abs`, `fp`, `ip` and such non analytic
    functions inside the integrand program. Other methods may behave
    better in this respect, but, on the other side, the double
    exponential method manages discontinuities (of the integrand or its
    derivatives) at the interval edges better and is usually way faster.

  As many other quadrature methods the double exponential algorithm does
    not manage well highly oscillatory integrands. Here, highly
    oscillatory means that the integrand changes sign many (or infinite)
    times along the integration interval.
  """

  tanhsinh = False            # True if tanh-sinh case
  expsinh = False             # True if exp-sinh case
  bpa2z = False               # True if limits in {(+/-inf, 0), (0, +/-inf), (+/-inf, -/+inf)}
  chg = False                 # True if limits = (+/-inf, b)

  if isnan(a) or isnan(b):
    return (NaN, NaN, 0, 0, 0, [])

  if a == b:
    return (0, 0, 0, 0, 0, [])

  if isfinite(a) and isfinite(b):
    # tanh-sinh case
    bpa2 = (b + a)/2          # centre point
    bma2 = (b - a)/2          # half interval
    tanhsinh = True
  elif isfinite(a) or isfinite(b):
    # exp-sinh case
    chg = isfinite(b)
    (bpa2, bma2) = (b, a) if chg else (a, b)
    bma2 = sign(bma2)
    bpa2z = almosteq(bpa2, 0.0)
    expsinh = True
  else:
    # sinh-sinh case
    bpa2 = 0
    bma2 = sign(b)
    bpa2z = True

  pi2 = pi()/2
  pi4 = pi()/4

  # convergence threshold
  eps = power(10, -mp.dps)
  thr = 10*sqrt(eps)
  if bpa2z:
    eps = power(10, -(mp.dps/2)**2)
  # maximum allowed level
  levelmax = int(round(log(mp.dps, 2)) + 1)   # + 2) also acceptable

  # maximum t
  if tanhsinh:
    tmax = 2*min(1, abs(bma2))
  elif bpa2z:
    tmax = sqrt(eps)
  else:
    tmax = abs(1/bpa2/2)
  tmax = ln(tmax/eps)
  if not tanhsinh:
    tmax *= 2
##  tmax = ln(tmax/pi2)
  exptmax = tmax/pi2

  s = 0                       # s is the computed integral
  h = 2                       # rectangle width
  tnfe = 0                    # Total Number of Function Evaluations
  q_lvl = []                  # computed value of integral at each level
  # progress thru levels
  for level in range(levelmax + 1):
    # sp = s at previous level
    sp = s*bma2*pi2*h
    if chg:
      sp = -sp
    if level:
      q_lvl.append(sp)
    wsl = 0                   # weigthed sum at this level
    h /= 2
    expt = exp(h)             # exp(t)
    exph = expt**2 if level else expt
    # walk abscissas
    while True:
      # try to avoid the computation of too many exp functions...
      iexpt = 1/expt
      cht = (expt + iexpt)/2
      # pi/2*sinh(t)
      pi2sh = pi4*(expt - iexpt)
      # weight and abscissa
      w = r = exp(pi2sh)
      if not expsinh:
        iexppi2sh = 1/r
        w += iexppi2sh        # 2cosh(pi2sh)
        r -= iexppi2sh        # 2sinh(pi2sh)
        r /= w if tanhsinh else 2
        w /= 2                # cosh(pi2sh)
      if tanhsinh:
        w = 1/w**2
      try:
        fpl = f(bpa2 + bma2*r)
      except ArithmeticError:
        fpl = 0
      p = fpl*w if isnormal(fpl) else 0
      try:
        fmi = f(bpa2 + bma2/r) if expsinh else f(bpa2 - bma2*r)
      except ArithmeticError:
        fmi = 0
      tnfe += 2
      p += (fmi/w if expsinh else fmi*w) if isnormal(fmi) else 0
      p *= cht
      wsl += p
      # next exp(t)
      expt *= exph
      # done with level?
      if expt > exptmax:
        break
      # early test (mainly for the sinh-sinh case)
      if abs(p) <= abs(eps*wsl):
        break
    # end of abscissa loop

    s += wsl
    # add the 1st series term
    if not level:
      try:
        s += f(bpa2 + bma2) if expsinh else f(bpa2)
      except ArithmeticError:
        pass
      tnfe += 1

    # converged?
    if not s or (abs(2*abs(wsl) - abs(s)) < abs(thr*s)):
      break
  # end of level loop

  # iteration done, apply constant coefficients
  s *= bma2*pi2*h
  if chg:
    s = -s
  q_lvl.append(s)

  # check for bad results
  err = abs(sp - s)
  if 10*err >= abs(s):
    err = abs(err) + abs(s)
    s = 0

  variant = 0 if tanhsinh else 1 if expsinh else 2

  return (s, err, tnfe, level, variant, q_lvl)


def main():
  if __name__ == '__main__':
    pass
