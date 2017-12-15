# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 21:04:26 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.exp(2 * x) - 1.


def L(x, i, xValues):
    xi = xValues[i]
    product = 1.
    for j, xj in enumerate(xValues):
        if j != i:
            product = product * (x - xj) / (xi - xj)
    return product


def Lagrange(x, xValues, fValues):
    val = 0.
    for i, fi in enumerate(fValues):
        val = val + L(x, i, xValues) * fi
    return val


xValues = np.array([1, 1.1, 1.2, 1.3, 1.4])
fValues = f(xValues)
p = Lagrange(1.25, xValues, fValues)
print('Function Value f(1.25)=', f(1.25))
print('Lagrange polyn p(1.25)=', p)

plt.close('all')
x = np.arange(1.0, 1.41, 0.01)
ps = Lagrange(x, xValues, fValues)
plt.plot(x, f(x), '.', label='f(x)')
plt.plot(xValues, fValues, '*')
plt.plot(x, ps, '-.', label='p(x)')
plt.legend()
plt.grid()
