# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import matplotlib.pyplot as plt
import interpolation


def f(x):
    return np.exp(-x)


def df(x):
    return -np.exp(-x)


def d2f(x):
    return f(x)


def err2(x):
    return np.exp(-x) * (x-0.25) *(x-0.5) *(x-0.75) /6


def err1(x):
    return np.exp(-x) * (0.5-0.25) *(0.5-0.75) /6

xValues = np.array([0.25, 0.5, 0.75])
fValues = f(xValues)
dValues = interpolation.differenceTable(fValues)


print('Actual df(0.5)=', df(0.5))
print('Numeric df(0.5)', interpolation.df(0.5, xValues, dValues))
print('Actual d2f(0.5)=', d2f(0.5))
print('Numeric d2f(0.5)', interpolation.d2f(0.5, xValues, dValues))

error = f(xValues[0]) / 6 * \
        np.abs((xValues[1]-xValues[0]) * (xValues[1] - xValues[2]))
print('Max error is', error)

x = np.arange(0.25, 0.76, 0.01)
plt.close('all')
#plt.plot(x, interpolation.df(x, xValues, dValues), label='arithmetic')
#plt.plot(x, df(x), label='actual')
plt.plot(x, err1(x), label='err1')
plt.plot(x, err2(x), label='err2')
plt.legend()
plt.grid()