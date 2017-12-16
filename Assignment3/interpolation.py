# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np


def df(x, xValues, dValues):
    h = xValues[1] - xValues[0]
    theta = (x - xValues[0]) / h
    return (dValues[1] + 0.5 * (2*theta - 1)*dValues[2]) / h


def d2f(x, xValues, dValues):
    h = xValues[1] - xValues[0]
    return dValues[2] / np.power(h, 2)


def DifferenceTable(fxValues):
    n = fxValues.size
    dValues = np.copy(fxValues)
    for i in range(1, n):
        for j in reversed(range(i, n)):
            dValues[j] = dValues[j] - dValues[j-1]
    return dValues


def dividedDifferenceTable(xValues, fxValues):
    n = fxValues.size
    dValues = np.copy(fxValues)
    for i in range(1, n):
        for j in reversed(range(i, n)):
            dValues[j] = (dValues[j] - dValues[j-1]) / \
                         (xValues[j] - xValues[j-i])
    return dValues


def dividedDifference(xValues, fxValues, s):
    val = 0.
    for i in range(s):
        xi = xValues[i]
        product = 1.
        for j, xj in enumerate(xValues):
            if j != i:
                product = product * (xi - xj)
        val = val + fxValues[i] / product
    return val


def NestedMultiplication(x, xValues, coeff):
    n = coeff.size
    y = coeff[n-1]
    for i in reversed(range(n - 1)):
        y = coeff[i] + (x - xValues[i]) * y
    return y