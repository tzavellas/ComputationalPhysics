# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 20:28:24 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import matplotlib.pyplot as plt
import interpolation

xValues = np.array([-1, 0., 1, 2, 3])
yValues = np.array([2., 1, 2, -7, 10])
coefficients = np.array([2., -1, 1, -2])

plt.close('all')
x = np.arange(-1, 3.1, 0.1)
y = interpolation.NestedMultiplication(x,
                                       xValues,
                                       coefficients)

plt.plot(xValues, yValues, 'ro')  # Plot points
plt.plot(x, y, '.', label='Newton 3rd')  # Plot 3rd degree Newton

newCoefficients = interpolation.addCoefficient(xValues,
                                               yValues,
                                               coefficients)
print('Order', newCoefficients.size - 1,
      'coefficent is', newCoefficients[-1])

y = interpolation.NestedMultiplication(x,
                                       xValues,
                                       newCoefficients)
plt.plot(x, y, '--', label='Newton 4th')  # Plot 4th degree Newton
plt.legend()
plt.grid()
