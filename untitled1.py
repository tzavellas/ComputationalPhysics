# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 14:10:12 2017

@author: extern.a.Tzavellas
"""

import numpy as np
import matplotlib.pyplot as plt


def fixedPoint(f, p0, tol, Nmax=None):
    i = 1
    p = p0
    while (True if (Nmax is None) else (i < Nmax)):
        p = f(p0)
        if (np.abs(p-p0) < tol):
            break
        i = i + 1
        p0 = p
    return {'p': p, 'nIter': i, 'found': True}


def g(x):
    return np.sqrt(-10 * np.cos(x))


def f(x):
    return np.power(x, 2) + 10 * np.cos(x)


p = -2
result = fixedPoint(g, p, 1e-4)
print('Solution = ', result['p'])
print('Iterations = ', result['nIter'])

x = np.arange(-5, 5, 0.1)
y = f(x)
plt.plot(x, y)
plt.grid()
