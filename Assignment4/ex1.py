# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 15:17:39 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import integration


def f(x):
    return np.power(x, 2) + np.power(x, 5)


IAnalytic = np.power(3, 3) / 3 + np.power(3, 6) / 6
print('Analytic')
print('Integral:', IAnalytic)

a = 0.
b = 3.
epsilon = 1e-2
trapezoid = integration.trapezoid(f, a, b, epsilon)
ITrapezoid = trapezoid['integral']
itTrapezoid = trapezoid['iterations']
print('\nTrapezoid')
print(trapezoid)
print('Deviation:', np.abs(IAnalytic - ITrapezoid))

simpson = integration.simpson(f, a, b, epsilon)
ISimpson = simpson['integral']
itSimpson = simpson['iterations']
print('\nSimpson 1/3')
print(simpson)
print('Deviation:', np.abs(IAnalytic - ISimpson))
