# Copyright (c) 2021, emece67 - MIT License


"""
A list of test cases to check the double exponential quadrature method.

The list is called `test_integral` and each item on it is a dict
  containing elements with keys:
    * 'a':  lower integration limit
    * 'b':  upper integration limit
    * 's':  true value of the integral
    * 'f':  a lambda expression taking one parameter `x` and returning
              the value of the function to be integrated
    * 'fs': a descriptive string

Can use mpmath and mp.dps must be set prior to importing this file.
"""


from mpmath import *


test_integral = [
  # Sinh-Sinh test functions  **************************************************
  {
    'a': +inf,
    'b': -inf,
    's': -2/sqrt(3),
    'f': lambda x: 1/(1 + x**2)/sqrt(3 + 3*x**2),
    'fs': "1/(1 + x**2)/sqrt(3 + 3*x**2)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()/12*ln(mpf('2/3')),
    'f': lambda x: x/(9*exp(x) + 4*exp(-x)),
    'fs': "x/(9*exp(x) + 4*exp(-x))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': ln(3)**2/4,
    'f': lambda x: x/(3 + exp(x))/(1 + exp(-x)),
    'fs': "x/(3 + exp(x))/(1 + exp(-x))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': ln(3)/3,
    'f': lambda x: x*exp(x)/(3 + exp(x))**2,
    'fs': "x*exp(x)/(3 + exp(x))**2"
  },
  {
    'a': -inf,
    'b': +inf,
    's': -ln(2),
    'f': lambda x: 1 - sqrt(2)*cosh(x)/sqrt(cosh(2*x)),
    'fs': "1 - sqrt(2)*cosh(x)/sqrt(cosh(2*x))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()/6*tan(pi()/3),
    'f': lambda x: exp(-x)*sinh(x)/sinh(3*x) if x != 0 else mpf('1/3'),
    'fs': "exp(-x)*sinh(x)/sinh(3*x) if x != 0 else mpf('1/3')"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()/3/sin(2*pi()/3),
    'f': lambda x: exp(-2*x)/(1 + exp(-3*x)),
    'fs': "exp(-2*x)/(1 + exp(-3*x))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': exp(mpf('1/9'))*sqrt(pi())/3,
    'f': lambda x: exp(-9*x**2 + 2*x),
    'fs': "exp(-9*x**2 + 2*x)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': (pi()/3/sin(2*pi()/3))**2,
    'f': lambda x: exp(2*x)*x/(exp(3*x) - 1) if x != 0 else mpf('1/3'),
    'fs': "exp(2*x)*x/(exp(3*x) - 1) if x != 0 else mpf('1/3')"
  },
  {
    'a': -inf,
    'b': +inf,
    's': mpf('3/2')*sqrt(pi()/2)*exp(mpf('9/2')),
    'f': lambda x: x*exp(-2*x**2 + 2*3*x),
    'fs': "x*exp(-2*x**2 + 2*3*x)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()**2/3,
    'f': lambda x: ln(abs((1 + 2*sqrt(1 + x**2))/(1 - 2*sqrt(1 + x**2))))/sqrt(1 + x**2),
    'fs': "ln(abs((1 + 2*sqrt(1 + x**2))/(1 - 2*sqrt(1 + x**2))))/sqrt(1 + x**2)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': mpf('2'),
    'f': lambda x: 1/cosh(x)**2,
    'fs': "1/cosh(x)**2"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi()/3)/6,
    'f': lambda x: x**2*exp(-3*x**2),
    'fs': "x**2*exp(-3*x**2)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi(),
    'f': lambda x: 1/(1 + x**2),
    'fs': "1/(1 + x**2)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': exp(mpf('4/9'))*sqrt(pi())/3,
    'f': lambda x: exp(-9*x**2 + 4*x),
    'fs': "exp(-9*x**2 + 4*x)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': (pi()/3/sin(2*pi()/3))**2,
    'f': lambda x: x*exp(2*x)/(exp(3*x) - 1) if x != 0 else mpf('1/3'),
    'fs': "x*exp(2*x)/(exp(3*x) - 1) if x != 0 else mpf('1/3')"
  },
  {
    'a': -inf,
    'b': +inf,
    's': ln(tan(3*pi()/8)/tan(2*pi()/8)),
    'f': lambda x: (exp(3*x) - exp(2*x))/x/(1 + exp(4*x)) if x != 0 else mpf('1/2'),
    'fs': "(exp(3*x) - exp(2*x))/x/(1 + exp(4*x)) if x != 0 else mpf('1/2')"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()/24*ln(mpf('4/3')),
    'f': lambda x: x/(9*exp(x) + 16*exp(-x)),
    'fs': "x/(9*exp(x) + 16*exp(-x))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': ln(3)**2/4,
    'f': lambda x: x/(3 + exp(x))/(1 + exp(-x)),
    'fs': "x/(3 + exp(x))/(1 + exp(-x))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': (pi()**2 + ln(3)**2)/8,
    'f': lambda x: x/(3 + exp(x))/(1 - exp(-x)) if x != 0 else mpf('1/4'),
    'fs': "x/(3 + exp(x))/(1 - exp(-x)) if x != 0 else mpf('1/4')"
  },
  {
    'a': -inf,
    'b': +inf,
    's': (pi()**2 + ln(3)**2)**2/16,
    'f': lambda x: x**3/(3 + exp(x))/(1 - exp(-x)) if x != 0 else 0,
    'fs': "x**3/(3 + exp(x))/(1 - exp(-x)) if x != 0 else 0"
  },
  {
    'a': -inf,
    'b': +inf,
    's': ln(3)/3,
    'f': lambda x: x*exp(x)/(3 + exp(x))**2,
    'fs': "x*exp(x)/(3 + exp(x))**2"
  },
  {
    'a': -inf,
    'b': +inf,
    's': ln(3)**2/2,
    'f': lambda x: x**2*(exp(x) - 3*exp(-x))/(3 + exp(x))**2/(1 + exp(-x))**2,
    'fs': "x**2*(exp(x) - 3*exp(-x))/(3 + exp(x))**2/(1 + exp(-x))**2"
  },
  {
    'a': -inf,
    'b': +inf,
    's': mpf('4/3')*sqrt(pi()/3)*exp(mpf('16/3')),
    'f': lambda x: x*exp(-3*x**2 + 8*x),
    'fs': "x*exp(-3*x**2 + 8*x)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': 1 - pi()/2/tan(pi()/2),
    'f': lambda x: sinh(x/2)**2/sinh(x)**2 if x!= 0 else mpf('1/4'),
    'fs': "sinh(x/2)**2/sinh(x)**2 if x!= 0 else mpf('1/4')"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()**2/3,
    'f': lambda x: x**2/sinh(x)**2 if x!= 0 else 1,
    'fs': "x**2/sinh(x)**2 if x!= 0 else 1"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()/10*tan(2*pi()/5),
    'f': lambda x: exp(-2*x)*sinh(2*x)/sinh(5*x) if x != 0 else mpf('4/10'),
    'fs': "exp(-2*x)*sinh(2*x)/sinh(5*x) if x != 0 else mpf('4/10')"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()*3/60,
    'f': lambda x: x*atan(3*x)/(x**2 + 9)/(x**2 + 9),
    'fs': "x*atan(3*x)/(x**2 + 9)/(x**2 + 9)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': 3*pi()/108,
    'f': lambda x: x*atan(3/x)/(x**2 + 9)/(x**2 + 9) if x != 0 else 0,
    'fs': "x*atan(3/x)/(x**2 + 9)/(x**2 + 9) if x != 0 else 0"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi())/3*exp(mpf('-1/9'))*sin(8),
    'f': lambda x: exp(-9*x**2)*sin(2*(x + 4)),
    'fs': "exp(-9*x**2)*sin(2*(x + 4))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi())/3*exp(mpf('-1/9'))*cos(8),
    'f': lambda x: exp(-9*x**2)*cos(2*(x + 4)),
    'fs': "exp(-9*x**2)*cos(2*(x + 4))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()*ln(1 + 6*abs(sin(5)) + 9),
    'f': lambda x: ln(9 - 2*3*x*cos(5) + x**2)/(1 + x**2),
    'fs': "ln(9 - 2*3*x*cos(5) + x**2)/(1 + x**2)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi()/3),
    'f': lambda x: exp(-3*x**2),
    'fs': "exp(-3*x**2)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi()/3)/12,
    'f': lambda x: x**4*exp(-3*x**2),
    'fs': "x**4*exp(-3*x**2)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi())/4*exp(mpf('-9/64'))*sin(15),
    'f': lambda x: exp(-16*x**2)*sin(3*(x + 5)),
    'fs': "exp(-16*x**2)*sin(3*(x + 5))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi())/4*exp(mpf('-9/64'))*cos(15),
    'f': lambda x: exp(-16*x**2)*cos(3*(x + 5)),
    'fs': "exp(-16*x**2)*cos(3*(x + 5))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi()*(cos(1) + sin(1))/4/exp(1),
    'f': lambda x: cos(x)/(x**4 + 4),
    'fs': "cos(x)/(x**4 + 4)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi()/3)*exp(mpf('-44/12')),
    'f': lambda x: exp(-(3*x**2 + 4*x + 5)),
    'fs': "exp(-(3*x**2 + 4*x + 5))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': 3*pi()/4,
    'f': lambda x: sin(x)**3/x**3 if x != 0 else 1,
    'fs': "sin(x)**3/x**3 if x != 0 else 1"
  },
  {
    'a': -inf,
    'b': +inf,
    's':  pi()*exp(-3*sin(3))*sin(12 + 3*cos(3)),
    'f': lambda x: (cos(9) - x*cos(12))/(1 - 2*x*cos(3) + x**2)*cos(3*x),
    'fs': "(cos(9) - x*cos(12))/(1 - 2*x*cos(3) + x**2)*cos(3*x)"
  },
  {
    'a': -inf,
    'b': +inf,
    's': sqrt(pi()),
    'f': lambda x: exp(-(x**2)),
    'fs': "exp(-(x**2))"
  },
  {
    'a': -inf,
    'b': +inf,
    's': pi(),
    'f': lambda x: x**2/(1 + 4*x + 3*x**2 - 4*x**3 - 2*x**4 + 2*x**5 + x**6),
    'fs': "x**2/(1 + 4*x + 3*x**2 - 4*x**3 - 2*x**4 + 2*x**5 + x**6)"
  },
  # Exp-Sinh test functions ****************************************************
  {
    'a': +inf,
    'b': mpf('-1'),
    's': -exp(mpf('5'))*sqrt(pi()/5),
    'f': lambda x: exp(-5*x)/sqrt(1 + x),
    'fs': "exp(-5*x)/sqrt(1 + x)"
  },
  {
    'a': mpf('1'),
    'b': +inf,
    's': pi()/cos(pi()/4),
    'f': lambda x: (x - 1)**(mpf('-1/4'))/x,
    'fs': "(x - 1)**(mpf('-1/4'))/x"
  },
  {
    'a': mpf('1'),
    'b': +inf,
    's': -pi()/3/sin(pi()/2)*sqrt(3),
    'f': lambda x: 1/(2 - 3*x)/(x - 1)**mpf('1/2'),
    'fs': "1/(2 - 3*x)/(x - 1)**mpf('1/2')"
  },
  {
    'a': mpf('3'),
    'b': +inf,
    's': pi(),
    'f': lambda x: 1/sqrt(x - 3)/(x - 2),
    'fs': "1/sqrt(x - 3)/(x - 2)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()/2,
    'f': lambda x: 1/(1 + x**2),
    'fs': "1/(1 + x**2)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': sqrt(pi()),
    'f': lambda x: exp(-x)/sqrt(x),
    'fs': "exp(-x)/sqrt(x)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': mpf('1/32'),
    'f': lambda x: x**2*exp(-4*x),
    'fs': "x**2*exp(-4*x)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': mpf('243/8'),
    'f': lambda x: (sqrt(x**2 + 9) - x)**3,
    'fs': "(sqrt(x**2 + 9) - x)**3"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': mpf('1/3'),
    'f': lambda x: exp(-3*x),
    'fs': "exp(-3*x)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': 1 - pi()**2/12,
    'f': lambda x: x*exp(-2*x)/(exp(-x) + 1),
    'fs': "x*exp(-2*x)/(exp(-x) + 1)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': ln(3)/3,
    'f': lambda x: ln(x)/(x + 3)**2,
    'fs': "ln(x)/(x + 3)**2"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()**2/12,
    'f': lambda x: ln(1 + exp(-x)),
    'fs': "ln(1 + exp(-x))"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': -ln(mpf('3/4'))/3,
    'f': lambda x: ln(1 + x)/(3*x + 4)**2,
    'fs': "ln(1 + x)/(3*x + 4)**2"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()**2/12,
    'f': lambda x: ln(1 + x**2)/x/(1 + x**2),
    'fs': "ln(1 + x**2)/x/(1 + x**2)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()*ln(2),
    'f': lambda x: ln((1 + x**2)/x)/(1 + x**2),
    'fs': "ln((1 + x**2)/x)/(1 + x**2)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': sqrt(pi()/3),
    'f': lambda x: exp(-3*x)/sqrt(x),
    'fs': "exp(-3*x)/sqrt(x)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi(),
    'f': lambda x: (ln(1 + 16*x**2) - ln(1 + 9*x**2))/x**2,
    'fs': "(ln(1 + 16*x**2) - ln(1 + 9*x**2))/x**2"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': sqrt(pi())/4,
    'f': lambda x: exp(-4*x**2),
    'fs': "exp(-4*x**2)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()**2/6 - 1,
    'f': lambda x: x*exp(-x)/(exp(x) - 1),
    'fs': "x*exp(-x)/(exp(x) - 1)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': ln(2) - mpf('1/2'),
    'f': lambda x: x/(1 + x**2)/sinh(pi()*x),
    'fs': "x/(1 + x**2)/sinh(pi()*x)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': -pi()/2,
    'f': lambda x: (1 - x**2)*ln(x)/(1 + x**2)**2,
    'fs': "(1 - x**2)*ln(x)/(1 + x**2)**2"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()**2/6,
    'f': lambda x: ln(1 + x)/x/(1 + x),
    'fs': "ln(1 + x)/x/(1 + x)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi(),
    'f': lambda x: ln(1 + x**2)/x**2,
    'fs': "ln(1 + x**2)/x**2"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': (mpf('1/2') + ln(2))*pi()/4,
    'f': lambda x: x*exp(-x)*sqrt(1 - exp(-2*x)),
    'fs': "x*exp(-x)*sqrt(1 - exp(-2*x))"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': 6,
    'f': lambda x: (sqrt(x**2 + 4) - x)**3,
    'fs': "(sqrt(x**2 + 4) - x)**3"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': mpf('3/32'),
    'f': lambda x: 1/(sqrt(x**2 + 4) + x)**3,
    'fs': "1/(sqrt(x**2 + 4) + x)**3"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': sqrt(pi()/2)/2*exp(-2*sqrt(6)),
    'f': lambda x: exp(-2*x**2 - 3/x**2),
    'fs': "exp(-2*x**2 - 3/x**2)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()**2/12 - mpf('3/4'),
    'f': lambda x: x*exp(-3*x)/(exp(-x) + 1),
    'fs': "x*exp(-3*x)/(exp(-x) + 1)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': 2*pi()**2/27,
    'f': lambda x: x*(1 - exp(-x))*exp(-x)/(exp(-3*x) + 1),
    'fs': "x*(1 - exp(-x))*exp(-x)/(exp(-3*x) + 1)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': 3*pi()/4*(ln(2) - mpf('7/12')),
    'f': lambda x: (x*exp(-2*x))/sqrt(exp(x) - 1),
    'fs': "(x*exp(-2*x))/sqrt(exp(x) - 1)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()/8*ln(12),
    'f': lambda x: ln(3*x)/(16 + x**2),
    'fs': "ln(3*x)/(16 + x**2)"
  },
  {
    'a': mpf('0'),
    'b': +inf,
    's': pi()**2/16,
    'f': lambda x: x*atan(x)/(1 + x**4),
    'fs': "x*atan(x)/(1 + x**4)"
  },
  # Tanh-Sinh test functions  **************************************************
  # easy functions
  {
    'a': mpf('1'),
    'b': mpf('0'),
    's': mpf('-1/4'),
    'f': lambda x: x*ln(1 + x),
    'fs': "x*ln(1 + x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': (pi() - 2 + 2*ln(mpf('2')))/12,
    'f': lambda x: x**2*atan(x),
    'fs': "x**2*atan(x)"
  },
  {
    'a': mpf('0'),
    'b': pi()/2,
    's': (exp(pi()/2) - 1)/2,
    'f': lambda x: exp(x)*cos(x),
    'fs': "exp(x)*cos(x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': 5*pi()**2/96,
    'f': lambda x: atan(sqrt(2 + x**2))/(1 + x**2)/sqrt(2 + x**2),
    'fs': "atan(sqrt(2 + x**2))/(1 + x**2)/sqrt(2 + x**2)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': pi(),
    'f': lambda x: 4/(x**2 + 1),
    'fs': "4/(x**2 + 1)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': 2/pi(),
    'f': lambda x: sin(pi()*x),
    'fs': "sin(pi()*x)"
  },
  {
    'a': mpf('0'),
    'b': pi(),
    's': pi()*asin(mpf('1/2')),
    'f': lambda x: mpf('1/2') if x == pi()/2 else ln(1 + cos(x)/2)/cos(x),
    'fs': "mpf('1/2') if x == pi()/2 else ln(1 + cos(x)/2)/cos(x)"
  },
  # Infinite derivative at one end
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('-4/9'),
    'f': lambda x: sqrt(x)*ln(x),
    'fs': "sqrt(x)*ln(x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('2/3'),
    'f': lambda x: sqrt(x),
    'fs': "sqrt(x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('2'),
    's': pi(),
    'f': lambda x: sqrt(x*(4 - x)),
    'fs': "sqrt(x*(4 - x))"
  },
  # Infinite derivative at both ends
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': pi()/4,
    'f': lambda x: sqrt(1 - x**2),
    'fs': "sqrt(1 - x**2)"
  },
  {
    'a': mpf('-1'),
    'b': mpf('1'),
    's': mpf('2.203345731824743771806893505443707840970711051705028816975706641754'),
    'f': lambda x: sqrt((1 - x**2)*(2 - x)),
    'fs': "sqrt((1 - x**2)*(2 - x))"
  },
  # Infinite derivative in the middle
  {
    'a': mpf('0'),
    'b': mpf('2'),
    's': mpf('4/3'),
    'f': lambda x: sqrt(abs(x - 1)),
    'fs': "sqrt(abs(x - 1))"
  },
  # Infinite derivative almost in the middle
  {
    'a': mpf('1/128'),
    'b': mpf('257/128'),
    's': ((129/mpf('128'))**mpf('3/2') + (127/mpf('128'))**mpf('3/2'))*2/3,
    'f': lambda x: sqrt(abs(x - 1)),
    'fs': "sqrt(abs(x - 1))"
  },
  # Infinite derivative at one end, infinite value at the other
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': 2*sqrt(pi())*gamma(mpf('3/4'))/gamma(mpf('1/4')),
    'f': lambda x: sqrt(x)/sqrt(1 - x**2),
    'fs': "sqrt(x)/sqrt(1 - x**2)"
  },
  {
    'a': mpf('0'),
    'b': pi()/2,
    's': pi()*sqrt(mpf('2'))/2,
    'f': lambda x: sqrt(tan(x)),
    'fs': "sqrt(tan(x))"
  },
  # Infinite value at one end
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('2'),
    'f': lambda x: ln(x)**2,
    'fs': "ln(x)**2"
  },
  {
    'a': mpf('0'),
    'b': pi()/2,
    's': -pi()*ln(mpf('2'))/2,
    'f': lambda x: ln(cos(x)),
    'fs': "ln(cos(x))"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('-1'),
    'f': lambda x: ln(x),
    'fs': "ln(x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('2'),
    'f': lambda x: 1/sqrt(x),
    'fs': "1/sqrt(x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': -si(mpf('1')),
    'f': lambda x: cos(x)*ln(x),
    'fs': "cos(x)*ln(x)"
  },
  # Infinite value at one end, but evaluates to 0
  {
    'a': mpf('0'),
    'b': e(),
    's': mpf('0'),
    'f': lambda x: ln(x),
    'fs': "ln(x)"
  },
  # Evaluates to 0
  {
    'a': mpf('-1'),
    'b': mpf('1'),
    's': mpf('0'),
    'f': lambda x: legendre(2, x),
    'fs': "legendre(2, x)"
  },
  {
    'a': mpf('-128'),
    'b': mpf('128'),
    's': mpf('0'),
    'f': lambda x: x**3*(x**2-47**2)*(x**2-88**2)*(x**2-117**2),
    'fs': "x**3*(x**2-47**2)*(x**2-88**2)*(x**2-117**2)"
  },
  {
    'a': -pi(),
    'b': pi(),
    's': mpf('0'),
    'f': lambda x: cos(x),
    'fs': "cos(x)"
  },
  # Bell shaped
  {
    'a': mpf('0'),
    'b': mpf('15'),
    's': pi()**4/15,
    'f': lambda x: x**3/(exp(x) - 1),
    'fs': "x**3/(exp(x) - 1)"
  },
  {
    'a': -10/sqrt(2),
    'b': 10/sqrt(2),
    's': sqrt(pi())*erf(10/sqrt(2)),
    'f': lambda x: exp(-x**2),
    'fs': "exp(-x**2)"
  },
  # Oscillatory
  {
    'a': -8*pi(),
    'b': 8.5*pi(),
    's': mpf('1'),
    'f': lambda x: cos(x),
    'fs': "cos(x)"
  },
  {
    'a': mpf('0'),
    'b': 2*pi(),
    's': mpf('-0.1976268077187172613672189604425462167176209935877442769219234953971'),
    'f': lambda x: ln(1 + x)*sin(10*x),
    'fs': "ln(1 + x)*sin(10*x)"
  },
  # Oscillatory, discontinuous
  {
    'a': mpf('0'),
    'b': mpf('6.4'),
    's': mpf('3.08'),
    'f': lambda x: frac(x),
    'fs': "frac(x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('7.4'),
    's': mpf('3.58'),
    'f': lambda x: frac(x),
    'fs': "frac(x)"
  },
  # Highly oscillatory, bounded
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('1/2'),
    'f': lambda x: cos(ln(x)),
    'fs': "cos(ln(x))"
  },
  {
    'a': mpf('0'),
    'b': mpf('2'),
    's': mpf('1.011239090533'),
    'f': lambda x: x**(pi()/4)*sin(pi()/(8 - 4*x)),
    'fs': "x**(pi()/4)*sin(pi()/(8 - 4*x))"
  },
  {
    'a': mpf('0'),
    'b': mpf('8'),
    's': si(mpf('1/8')) + 8*cos(mpf('1/8')) - pi()/2,
    'f': lambda x: cos(1/(8 - x)),
    'fs': "cos(1/(8 - x))"
  },
  # Highly oscillatory, unbounded
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': pi()/2 - si(1),
    'f': lambda x: sin(1/x)/x,
    'fs': "sin(1/x)/x"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('0.3233674314850376872853210319016270981575283855468855736552649593789637419574556193621185898284113307'),
    'f': lambda x: cos(ln(x)/x)/x,
    'fs': "cos(ln(x)/x)/x"
  },
  # Kahan's pathological
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': 2 + digamma(1) - ln(4),
    'f': lambda x: 2*x**2/((x-1)*(x+1)) - x/ln(x),
    'fs': "2*x**2/((x-1)*(x+1)) - x/ln(x)"
  },
  # A trap
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': +inf,
    'f': lambda x: 1/x,
    'fs': "1/x"
  },
  # BI nightmare poly
  {
    'a': mpf('-128'),
    'b': mpf('128'),
    's': mpf('8588888219586658304/315'),
    'f': lambda x: x**2*(x**2 - 47**2)*(x**2 - 88**2)*(x**2 - 117**2),
    'fs': "x**2*(x**2 - 47**2)*(x**2 - 88**2)*(x**2 - 117**2)"
  },
  {
    'a': mpf('-128'),
    'b': mpf('128'),
    's': mpf('11816497947871281152/315'),
    'f': lambda x: x**2*(x**2 - 46**2)*(x**2 - 87**2)*(x**2 - 116**2),
    'fs': "x**2*(x**2 - 46**2)*(x**2 - 87**2)*(x**2 - 116**2)"
  },
  # these are corner cases  ****************************************************
  {
    'a': mpf('-1'),
    'b': mpf('1'),
    's': mpf('2/8'),
    'f': lambda x: abs(x)*x**6,
    'fs': "abs(x)*x**6"
  },
  {
    'a': -inf,
    'b': inf,
    's': pi(),
    'f': lambda x: sin(x)/x if x != 0 else 1,
    'fs': "sin(x)/x if x != 0 else 1"
  },
  {
    'a': 0,
    'b': inf,
    's': pi()/2,
    'f': lambda x: sin(x)/x if x != 0 else 1,
    'fs': "sin(x)/x if x != 0 else 1"
  },
  {
    'a': +inf,
    'b': mpf('0'),
    's': -pi()/2,
    'f': lambda x: 1/(1 + x**2),
    'fs': "1/(1 + x**2)"
  },
  {
    'a': -inf,
    'b': mpf('1'),
    's': exp(mpf('5'))*sqrt(pi()/5),
    'f': lambda x: exp(5*x)/sqrt(1 - x),
    'fs': "exp(5*x)/sqrt(1 - x)"
  },
  {
    'a': mpf('1'),
    'b': -inf,
    's': -exp(mpf('5'))*sqrt(pi()/5),
    'f': lambda x: exp(5*x)/sqrt(1 - x),
    'fs': "exp(5*x)/sqrt(1 - x)"
  },
  {
    'a': -inf,
    'b': mpf('0'),
    's': -exp(mpf('5'))*sqrt(pi())*(erf(sqrt(mpf('5'))) - 1)/sqrt(mpf('5')),
    'f': lambda x: exp(5*x)/sqrt(1 - x),
    'fs': "exp(5*x)/sqrt(1 - x)"
  },
  {
    'a': mpf('0'),
    'b': -inf,
    's': exp(mpf('5'))*sqrt(pi())*(erf(sqrt(mpf('5'))) - 1)/sqrt(mpf('5')),
    'f': lambda x: exp(5*x)/sqrt(1 - x),
    'fs': "exp(5*x)/sqrt(1 - x)"
  },
  {
    'a': inf,
    'b': mpf('-1'),
    's': -exp(mpf('5'))*sqrt(pi()/5),
    'f': lambda x: exp(-5*x)/sqrt(1 + x),
    'fs': "exp(-5*x)/sqrt(1 + x)"
  },
  {
    'a': mpf('-1'),
    'b': inf,
    's': exp(mpf('5'))*sqrt(pi()/5),
    'f': lambda x: exp(-5*x)/sqrt(1 + x),
    'fs': "exp(-5*x)/sqrt(1 + x)"
  },
  {
    'a': inf,
    'b': mpf('0'),
    's': exp(mpf('5'))*sqrt(pi())*(erf(sqrt(mpf('5'))) - 1)/sqrt(mpf('5')),
    'f': lambda x: exp(-5*x)/sqrt(1 + x),
    'fs': "exp(-5*x)/sqrt(1 + x)"
  },
  {
    'a': mpf('0'),
    'b': inf,
    's': -exp(mpf('5'))*sqrt(pi())*(erf(sqrt(mpf('5'))) - 1)/sqrt(mpf('5')),
    'f': lambda x: exp(-5*x)/sqrt(1 + x),
    'fs': "exp(-5*x)/sqrt(1 + x)"
  },
  {
    'a': -inf,
    'b': -1,
    's': pi()/4,
    'f': lambda x: 1/(x**2 + 1),
    'fs': "1/(x**2 + 1)"
  },
  {
    'a': -inf,
    'b': inf,
    's': 0,
    'f': lambda x: (x - 1)/((x - 1)**2 + 1)**2,
    'fs': "(x - 1)/((x - 1)**2 + 1)**2"
  },
  {
    'a': -inf,
    'b': 0,
    's': mpf('-1/4'),
    'f': lambda x: (x - 1)/((x - 1)**2 + 1)**2,
    'fs': "(x - 1)/((x - 1)**2 + 1)**2"
  },
  {
    'a': -inf,
    'b': 1,
    's': mpf('-1/2'),
    'f': lambda x: (x - 1)/((x - 1)**2 + 1)**2,
    'fs': "(x - 1)/((x - 1)**2 + 1)**2"
  },
  {
    'a': -inf,
    'b': 2,
    's': mpf('-1/4'),
    'f': lambda x: (x - 1)/((x - 1)**2 + 1)**2,
    'fs': "(x - 1)/((x - 1)**2 + 1)**2"
  },
  {
    'a': 0,
    'b': 2,
    's': 0,
    'f': lambda x: (x - 1)/((x - 1)**2 + 1)**2,
    'fs': "(x - 1)/((x - 1)**2 + 1)**2"
  },
  {
    'a': 0,
    'b': 1,
    's': mpf('-1/4'),
    'f': lambda x: (x - 1)/((x - 1)**2 + 1)**2,
    'fs': "(x - 1)/((x - 1)**2 + 1)**2"
  },
  {
    'a': +inf,
    'b': -inf,
    's': -mpf('2e-30')/sqrt(3),
    'f': lambda x: mpf('1e-30')/(1 + x**2)/sqrt(3 + 3*x**2),
    'fs': "mpf('1e-30')/(1 + x**2)/sqrt(3 + 3*x**2)"
  },
  {
    'a': 0,
    'b': pi(),
    's': mpf('2e-30'),
    'f': lambda x: mpf('1e-30')*sin(x),
    'fs': "mpf('1e-30')*sin(x)"
  },
  {
    'a': +inf,
    'b': mpf('-1'),
    's': -mpf('1e-30')*exp(mpf('5'))*sqrt(pi()/5),
    'f': lambda x: mpf('1e-30')*exp(-5*x)/sqrt(1 + x),
    'fs': "mpf('1e-30')*exp(-5*x)/sqrt(1 + x)"
  },
  {
    'a': pi()*4,
    'b': pi()*5,
    's': mpf('0.002134240777991131240110232369370650561940890093246454073400113711660569572398558456002767541523306133688281981489760809834964047628791269'),
    'f': lambda x: sin(x)**exp(x),
    'fs': "sin(x)**exp(x)"
  },
  {
    'a': pi()*4.25,
    'b': pi()*5,
    's': mpf('0.002134240777991131240110232369370650561940890093246454073400113711660569572398558456002767541523306133688281981489760809834964047628791269'),
    'f': lambda x: sin(x)**exp(x),
    'fs': "sin(x)**exp(x)"
  },
  {
    'a': -pi(),
    'b': pi(),
    's': mpf('0'),
    'f': lambda x: mpf('0'),
    'fs': "mpf('0')"
  },
  # robve **********************************************************************
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('1/12'),
    'f': lambda x: x**3 - 2*x**2 + x,
    'fs': "x**3 - 2*x**2 + x"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': ln(2),
    'f': lambda x: 1/(1 + x),
    'fs': "1/(1 + x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': pi(),
    'f': lambda x: 4/(1 + x*x),
    'fs': "4/(1 + x*x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('1'),
    'f': lambda x: acos(x),
    'fs': "acos(x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('0.9460830703671830149413533138231796578123379547381117904714547735666870365407979180887021330817407112'),
    'f': lambda x: sin(x)/x,
    'fs': "sin(x)/x"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('1.198140234735592207439922492280323878227212663215651558263674952946405214143915670835885556489793389'),
    'f': lambda x: sqrt(x/(1 - x**2)),
    'fs': "sqrt(x/(1 - x**2))"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('2'),
    'f': lambda x: ln(x)**2,
    'fs': "ln(x)**2"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('2'),
    'f': lambda x: 1/sqrt(x),
    'fs': "1/sqrt(x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('2'),
    'f': lambda x: 1/sqrt(1 - x),
    'fs': "1/sqrt(1 - x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('5'),
    'f': lambda x: x**mpf('-0.8'),
    'fs': "x**mpf('-0.8')"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('5'),
    'f': lambda x: (1 - x)**mpf('-0.8'),
    'fs': "(1 - x)**mpf('-0.8')"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('1.669253683348146372562859465598093617987986026980694004899654740207363985419052823739382320702550648'),
    'f': lambda x: 1/sqrt(sin(pi()*x)),
    'fs': "1/sqrt(sin(pi()*x))"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('3.604250526330089151536169815574269902046879077863075298569203952881902348266779199510575127257374918'),
    'f': lambda x: sin(pi()*x)**mpf('-0.8'),
    'fs': "sin(pi()*x)**mpf('-0.8')"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': sqrt(pi()),
    'f': lambda x: 1/sqrt(-ln(x)),
    'fs': "1/sqrt(-ln(x))"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': sqrt(pi()),
    'f': lambda x: 1/sqrt(-ln(1 - x)),
    'fs': "1/sqrt(-ln(1 - x))"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('0'),
    'f': lambda x: sin(pi()*x*40),
    'fs': "sin(pi()*x*40)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': atan(5)/5,
    'f': lambda x: 1/(1 + 25*x**2),
    'fs': "1/(1 + 25*x**2)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': 5*atan(mpf('1/5')),
    'f': lambda x: 1/(1 + mpf('0.04')*x**2),
    'fs': "1/(1 + mpf('0.04')*x**2)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': sqrt(2)/3,
    'f': lambda x: sqrt(abs(x - mpf('0.5'))),
    'fs': "sqrt(abs(x - mpf('0.5')))"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('4.5'),
    'f': lambda x: floor(10*x),
    'fs': "floor(10*x)"
  },
  {
    'a': mpf('0'),
    'b': mpf('1'),
    's': mpf('0.5'),
    'f': lambda x: 10*x - floor(10*x),
    'fs': "10*x - floor(10*x)"
  }
]


def main():
  if __name__ == '__main__':
    pass
