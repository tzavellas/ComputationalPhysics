# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np


def dividedDifferenceTable(xValues, fxValues):
    n = xValues.size
    d = np.copy(fxValues)
    for i in range(1, n):
        for j in reversed(range(i, n)):
            d[j] = (d[j] - d[j-1]) / (xValues[j] - xValues[j-i])
    return d


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