# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 20:28:24 2017

@author: tzave
"""

import numpy as np
import matplotlib.pyplot as plt
import interpolation


def f(x):
    return x + np.sin(x)


def er(x):
    return x * (x - .5) * (x - 1) * np.sin(x) / np.math.factorial(4)


xValues = np.array([0, 0.5, 1, 1.5])
yValues = f(xValues)

coefficients = interpolation.dividedDifferenceTable(xValues, yValues)
p = interpolation.NestedMultiplication(0, xValues, coefficients)

plt.close('all')
x = np.arange(0, 1.6, 0.1)
y = interpolation.NestedMultiplication(x, xValues, coefficients)
plt.figure(1)
plt.plot(x, f(x), '.', label='x + sinx')
plt.plot(x, y, '--', label='Newton')
plt.plot(xValues, yValues, '*')
plt.legend()
plt.grid()


x = np.arange(0, np.pi/2, 0.01)
err = np.abs(er(x))
plt.figure(2)
plt.plot(x, err, label='|f(x)-p(x)|')
plt.legend()
plt.grid()
error = (np.power(np.pi, 3) - 3 * np.power(np.pi, 2) + 2 * np.pi) / 192
print('Theoretical Error:', error)
print('Graphical Error:', max(err))
