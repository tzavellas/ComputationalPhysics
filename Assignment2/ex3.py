# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 20:28:24 2017

@author: tzave
"""

import numpy as np
import matplotlib.pyplot as plt


def NewtonRaphson(x0, f, df, tol, N=None):
    i = 1
    p0 = x0
    while (True if N is None else (i < N)):
        p = p0 - f(p0)/df(p0)
        if np.abs(p-p0) < tol:
            break
        i = i + 1
        p0 = p
    return {'p': p, 'iterations': i}


def RegulaFalsi(a, b, f, tol, N=None):
    i = 1
    p0 = a
    p1 = b
    q0 = f(a)
    q1 = f(b)
    while (True if N is None else (i < N)):
        p = p1 - q1 * (p1 - p0)/(q1 - q0)
        if np.abs(p - p1) < tol:
            break
        i = i + 1
        p0 = p1
        q0 = q1
        p1 = p
        q1 = f(p)
    return {'p': p, 'iterations': i}


def Bisection(a, b, f, tol, N=None):
    i = 1
    c_old = 0
    while (True if N is None else (i < N)):
        c = (a + b) / 2
        if np.isclose(f(c), 0) or (i > 1 and (np.abs(c-c_old) < tol)):
            break
        i = i + 1
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
        c_old = c
    return {'p': c, 'iterations': i}


def f(x):
    return np.power(x, 4) + 2 * np.power(x, 3) - 5 * np.power(x, 2) - 7 * x - 5


def df(x):
    return 4 * np.power(x, 3) + 6 * np.power(x, 2) - 10 * x - 7


x0 = 1.5
roots = np.zeros(3)

result = NewtonRaphson(x0, f, df, 1e-5)
roots[0] = result['p']
iterations = result['iterations']
print('Newton Raphson Method: Root ', roots[0], ' found after ', iterations, ' iterations')

result = RegulaFalsi(1.5, 3, f, 1e-5)
roots[1] = result['p']
iterations = result['iterations']
print('RegulaFalsi Method: Root ', roots[1], ' found after ', iterations, ' iterations')

result = Bisection(1.5, 3, f, 1e-5)
roots[2] = result['p']
iterations = result['iterations']
print('Bisection Method: Root ', roots[2], ' found after ', iterations, ' iterations')

x = np.arange(1.5, 3, 0.1)
plt.close('all')
plt.figure(1)
plt.plot(x, f(x), roots[0], f(roots[0]), '*r', roots[1], f(roots[1]), '*b', roots[2], f(roots[2]), '*k')
plt.grid()

x = np.arange(1.5, 3, 0.01)
plt.figure(2)
plt.plot(x, df(x), label = 'df')
plt.legend()
plt.grid()
