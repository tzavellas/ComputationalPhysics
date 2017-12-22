# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 14:51:41 2017

@author: Anastasios Tzavellas
"""


import numpy as np


def trapezoid(f, a, b, epsilon=1e-3):
    n = 1
    err = 1.
    integralOld = 0.
    while(err > epsilon):
        h = (b - a)/n
        sumf = 0.
        if n != 1:
            xValues = np.arange(a, b, h)
            fValues = f(xValues)
            for i in range(1, n):
                sumf = sumf + fValues[i]
        integralNew = h/2 * (f(a) + 2*sumf + f(b))
        if n != 1:
            err = np.abs(integralNew - integralOld)
        integralOld = integralNew
        n = 2 * n
    return {'integral': integralNew, 'iterations': np.log2(n).astype(int)}


def simpson(f, a, b, epsilon=1e-3):
    n = 2
    err = 1.
    integralOld = 0.
    while(err > epsilon):
        h = (b - a)/n
        sumfEven = 0.
        sumfOdd = 0.
        xValues = np.arange(a, b, h)
        fValues = f(xValues)
        for i in range(1, n):
            if i % 2 == 0:
                sumfEven = sumfEven + fValues[i]
            else:
                sumfOdd = sumfOdd + fValues[i]
        integralNew = h/3 * (f(a) + 2*sumfEven + 4*sumfOdd + f(b))
        if n != 2:
            err = np.abs(integralNew - integralOld)
        integralOld = integralNew
        n = 2 * n
    return {'integral': integralNew, 'iterations': np.log2(n).astype(int)}
