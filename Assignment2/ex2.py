# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 20:28:24 2017

@author: tzave
"""

import numpy as np
import matplotlib.pyplot as plt


def FixedPoint(x0, g, tol, N=None):
    i = 1
    p0 = x0
    while (True if N is None else (i < N)):
        p = g(p0)
        if np.abs(p-p0) < tol:
            break
        i = i + 1
        p0 = p
    return {'p': p, 'iterations': i}


def f(x):
    return np.power(x, 3) + -x - 1


def g1(x):
    return np.power(x + 1, 1/3)


def dg1(x):
    return np.power(x + 1, -2/3) / 3


x0 = 1
result = FixedPoint(x0, g1, 1e-5)
root = result['p']
iterations = result['iterations']
print('Root found after ', iterations, ' iterations')

x = np.arange(1, 2, 0.1)
plt.close('all')
plt.figure(1)
plt.plot(x, f(x), root, f(root), '*r')
plt.grid()

x = np.arange(1, 2, 0.01)
plt.figure(2)
plt.plot(x, dg1(x), label = 'dg1')
plt.plot(x, np.ones(x.shape), label = 'y = 1')
plt.plot(x, -np.ones(x.shape), label = 'y = -1')
plt.ylim(-1.5, 1.5)
plt.legend()
plt.grid()