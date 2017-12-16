# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import interpolation


def f(x):
    return np.exp(-x)


def df(x):
    return -np.exp(-x)


def d2f(x):
    return f(x)


xValues = np.array([0.25, 0.5, 0.75])
fValues = f(xValues)
dValues = interpolation.DifferenceTable(fValues)


print('Actual df(0.5)=', df(0.5))
print('Numeric df(0.5)', interpolation.df(0.5, xValues, dValues))
print('Actual d2f(0.5)=', d2f(0.5))
print('Numeric d2f(0.5)', interpolation.d2f(0.5, xValues, dValues))

error = f(xValues[0]) / 6 * \
        np.abs((xValues[1]-xValues[0]) * (xValues[1] - xValues[2]))
print('Max error is', error)
