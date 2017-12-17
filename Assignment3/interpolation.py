# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import scipy.integrate as integrate


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


def lSquaresRightHand(f, order, a, b):
    n = order + 1
    y = np.zeros(n)
    for i in range(n):
        y[i] = integrate.quad(lambda x: f(x) * np.power(x, i), a, b)[0]
        y[i] = 0. if np.isclose(y[i], 0., atol=1e-8) else y[i]
    return y


def VandermondeMatrix(order, a, b):
    n = order + 1
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = integrate.quad(lambda x: np.power(x, i+j), a, b)[0]
    return A
