# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 15:08:41 2017

@author: Anastasios Tzavellas
"""

import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import diffEquation as de


def f(t, y):
    g = 10
    c = g/4
    A = np.array([[0, 1, 0, 0],
                  [0, -c, 0, 0],
                  [0, 0, 0, 1],
                  [0, 0, -c, 0]])
    b = np.array([0, 0, 0, -g])
    return np.dot(A, y) + b


def RangeFun(phi, V0, h, f, solver):
    Vx = V0 * np.cos(phi)
    Vy = V0 * np.sin(phi)
    Y0 = np.array([0., Vx, 0., Vy])
    solution = de.solve(Y0, h, f, solver)
    x = solution['x']
    return x[x.size - 1]


plt.close('all')
V0 = 10.  # Initial velocity
h = 0.01

# Find maximum using Brent's Algorithm
g = minimize_scalar(lambda phi: -RangeFun(phi,
                                          V0,
                                          h,
                                          f,
                                          de.rungeKutta4),
                    bracket=(0, np.pi/2))
print("Max Range @%.2f deg: %.4f" % (np.rad2deg(g['x']),
                                     -g['fun']))
plt.figure(1)
phi_ = np.arange(0.,
                 np.pi/2,
                 0.0091)  # range of angles to find range
Ranges = list()
i = 0
for phi in phi_:
    Ranges.insert(i, RangeFun(phi,
                              V0,
                              h,
                              f,
                              de.rungeKutta4))
    i = i + 1
plt.plot(np.rad2deg(phi_),
         np.array(Ranges))  # Plot range vs angle
plt.ylabel('Range [m]')
plt.xlabel('Angle [deg]')
plt.grid()

# Exhaustive Search - Assumes range functions has one maximum
Range = 0.
Ranges = list()
Y = list()
for i, phi in enumerate(phi_):
    Vx = V0 * np.cos(phi)
    Vy = V0 * np.sin(phi)
    Y0 = np.array([0., Vx, 0., Vy])
    solution = de.solve(Y0,
                        h,
                        f,
                        de.rungeKutta4)
    Y.insert(i, solution)
    x = solution['x']
    y = solution['y']
    Ranges.insert(i, x[x.size - 1])

maxRange = max(Ranges)  # Find max element of calculated values
i = Ranges.index(maxRange)  # Find phi corresponding to max element
maxPhi = np.rad2deg(phi_[i])

print("Max Range @%.2f deg: %.4f" % (maxPhi,
                                  maxRange))
plt.figure(2)
plt.plot(Y[i]['x'],
         Y[i]['y'],
         label="%.2f deg" % maxPhi)
plt.ylabel('y [m]')
plt.xlabel('x [m]')
plt.grid()
plt.legend()
