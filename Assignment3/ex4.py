# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt


def f1(x):
    return np.power(x, 2) - 2 * x + 3


def f2(x):
    return np.cos(np.pi*x)


def f3(x):
    return np.exp(-x)


def lSquaresRightHand(f, order, a, b):
    n = order + 1
    y = np.zeros(n)
    for i in range(n):
        y[i] = integrate.quad(lambda x: f(x) * np.power(x, i), a, b)[0]
        y[i] = 0. if np.isclose(y[i], 0., atol=1e-8) else y[i]
    return y


def lSquaresMatrix(order, a, b):
    n = order + 1
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = integrate.quad(lambda x: np.power(x, i+j), a, b)[0]
    return A


a = 0
b = 1
order = 2

x = np.arange(a, b, 0.01)
A = lSquaresMatrix(order, a, b)

plt.close('all')

B1 = lSquaresRightHand(f1, order, a, b)
y1 = np.linalg.solve(A, B1)
p1 = np.poly1d(np.flip(y1, 0))  # Least Squares polynomial
plt.figure(1)
plt.plot(x, f1(x), label='x^2 -2x +3')
plt.plot(x, p1(x), label='polynomial')
plt.legend()
plt.grid()


B2 = lSquaresRightHand(f2, order, a, b)
y2 = np.linalg.solve(A, B2)
p2 = np.poly1d(np.flip(y2, 0))  # Least Squares polynomial
plt.figure(2)
plt.plot(x, f2(x), label='cos(pi x)')
plt.plot(x, p2(x), label='polynomial')
plt.legend()
plt.grid()

B3 = lSquaresRightHand(f3, order, a, b)
y3 = np.linalg.solve(A, B3)
p3 = np.poly1d(np.flip(y3, 0))  # Least Squares polynomial
plt.figure(3)
plt.plot(x, f3(x), label='exp(-x)')
plt.plot(x, p3(x), label='polynomial')
plt.legend()
plt.grid()
