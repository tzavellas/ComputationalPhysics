# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import matplotlib.pyplot as plt
import interpolation


def f1(x):
    """Evaluates x^2-2x+3"""
    return np.polyval([3, -2, 1], x)


def f2(x):
    """Evaluates cos(pi * x)"""
    return np.cos(np.pi*x)


def f3(x):
    """Evaluates e^(-x)"""
    return np.exp(-x)


def Legendre(x, n):
    """Evaluates Legendre polynomial
    of n degree"""
    c = np.zeros(n + 1)
    c[n] = 1.
    return np.polynomial.legendre.legval(x, c)


def wLegendre(x):
    """Weight function for Legendre
    polynomials used in Least Squares"""
    return np.ones(x)


a = 0
b = 1
order = 2

x = np.arange(a, b, 0.01)

plt.close('all')

p1 = interpolation.lSquaresMonomial(f1, order, a, b)
a1 = interpolation.lSquaresOrthogonal(f1,
                                      wLegendre,
                                      Legendre,
                                      3,
                                      -1,
                                      1)  # Base Legendre polynomials
print("Monomial:", p1)
print("Legendre:", a1)
plt.figure(1)
plt.plot(x, f1(x), label='x^2 -2x +3')
plt.plot(x, p1(x), '--', label='monomial')
plt.plot(x, a1(x), '.', label='Legendre')
plt.legend()
plt.grid()

p2 = interpolation.lSquaresMonomial(f2, order, a, b)
a2 = interpolation.lSquaresOrthogonal(f2,
                                      wLegendre,
                                      Legendre,
                                      3,
                                      -1,
                                      1)  # Base Legendre polynomials
print("\nMonomial:", p2)
print("Legendre:", a2)
plt.figure(2)
plt.plot(x, f2(x), label='cos(pi x)')
plt.plot(x, p2(x), '--', label='monomial')
plt.plot(x, a2(x), '.', label='Legendre')
plt.legend()
plt.grid()

p3 = interpolation.lSquaresMonomial(f3, order, a, b)
a3 = interpolation.lSquaresOrthogonal(f3,
                                      wLegendre,
                                      Legendre,
                                      3,
                                      -1,
                                      1)  # Base Legendre polynomials
print("\nMonomial:", p3)
print("Legendre:", a3)
plt.figure(3)
plt.plot(x, f3(x), label='exp(-x)')
plt.plot(x, p3(x), '--', label='monomial')
plt.plot(x, a3(x), '.', label='Legendre')
plt.legend()
plt.grid()
