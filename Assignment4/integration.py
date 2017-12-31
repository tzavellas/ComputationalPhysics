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
    fa = f(a)
    fb = f(b)
    evaluations = 2
    while(err > epsilon):
        h = (b - a)/n
        sumf = 0.
        if n != 1:
            xValues = np.arange(a, b, h)
            fValues = f(xValues)
            evaluations = evaluations + fValues.size
            for i in range(1, n):
                sumf = sumf + fValues[i]
        integralNew = h/2 * (fa + 2*sumf + fb)
        if n != 1:
            err = np.abs(integralNew - integralOld)
        integralOld = integralNew
        n = 2 * n
    return {'integral': integralNew,
            'iterations': np.log2(n).astype(int),
            'evaluations': evaluations}


def simpson(f, a, b, epsilon=1e-3):
    n = 2
    err = 1.
    integralOld = 0.
    evaluations = 2
    fa = f(a)
    fb = f(b)
    while(err > epsilon):
        h = (b - a)/n
        sumfEven = 0.
        sumfOdd = 0.
        xValues = np.arange(a, b, h)
        fValues = f(xValues)
        evaluations = evaluations + fValues.size
        for i in range(1, n):
            if i % 2 == 0:
                sumfEven = sumfEven + fValues[i]
            else:
                sumfOdd = sumfOdd + fValues[i]
        integralNew = h/3 * (fa + 2*sumfEven + 4*sumfOdd + fb)
        if n != 2:
            err = np.abs(integralNew - integralOld)
        integralOld = integralNew
        n = 2 * n
    return {'integral': integralNew,
            'iterations': np.log2(n).astype(int),
            'evaluations': evaluations}
