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
A = interpolation.normalEquations(order, a, b)

plt.close('all')

B1 = interpolation.lSquaresRightHand(f1, order, a, b)
y1 = np.linalg.solve(A, B1)  # Solve with base {1, x, x^2, ...}
p1 = np.poly1d(np.flip(y1, 0))  # LS polynomial
a1 = interpolation.lSquaresOrthogonal(f1,
                                      wLegendre,
                                      Legendre,
                                      3,
                                      -1,
                                      1)  # Base Legendre polynomials
plt.figure(1)
plt.plot(x, f1(x), label='x^2 -2x +3')
plt.plot(x, p1(x), '--', label='monomial')
plt.plot(x, np.polynomial.legendre.legval(x, a1), '.', label='Legendre')
plt.legend()
plt.grid()


B2 = interpolation.lSquaresRightHand(f2, order, a, b)
y2 = np.linalg.solve(A, B2)
p2 = np.poly1d(np.flip(y2, 0))  # Least Squares polynomial
a2 = interpolation.lSquaresOrthogonal(f2,
                                      wLegendre,
                                      Legendre,
                                      3,
                                      -1,
                                      1)  # Base Legendre polynomials
plt.figure(2)
plt.plot(x, f2(x), label='cos(pi x)')
plt.plot(x, p2(x), '--', label='monomial')
plt.plot(x, np.polynomial.legendre.legval(x, a2), '.', label='Legendre')
plt.legend()
plt.grid()

B3 = interpolation.lSquaresRightHand(f3, order, a, b)
y3 = np.linalg.solve(A, B3)
p3 = np.poly1d(np.flip(y3, 0))  # Least Squares polynomial
a3 = interpolation.lSquaresOrthogonal(f3,
                                      wLegendre,
                                      Legendre,
                                      3,
                                      -1,
                                      1)  # Base Legendre polynomials
plt.figure(3)
plt.plot(x, f3(x), label='exp(-x)')
plt.plot(x, p3(x), '--', label='monomial')
plt.plot(x, np.polynomial.legendre.legval(x, a3), '.', label='Legendre')
plt.legend()
plt.grid()
