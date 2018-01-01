# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 17:32:05 2017

@author: tzave
"""
import numpy as np


def eulerStep(ti, yi, f, h=0.1):
    return yi + h * f(ti, yi)


def rungeKutta4(ti, yi, f, h=0.1):
    k1 = f(ti, yi)
    k2 = f(ti + h/2, yi + h/2 * k1)
    k3 = f(ti + h/2, yi + h/2 * k2)
    k4 = f(ti + h, yi + h * k3)
    return yi + h/6 * (k1 + 2*k2 + 2*k3 + k4)


def solve(Y0, f, method, h=0.1, tmin=0, tmax=None, N=None):
    t = np.array([])  # store time
    x = np.array([])  # store x coord
    y = np.array([])  # store y coord
    vx = np.array([])  # store x velocity
    vy = np.array([])  # store y velocity
    i = 0
    Y = Y0
    while 1:
        x = np.insert(x, i, Y[0])  # store all previous step values
        y = np.insert(y, i, Y[2])
        vx = np.insert(vx, i, Y[1])
        vy = np.insert(vy, i, Y[3])
        t = np.insert(t, i, i * h)
        Y = method(t[i], Y, f, h)  # Next step
        if (Y[2] < 0.) or (i == N):  # break if max iterations
            break  # or if y<.0
        i = i + 1
    return {'x': x,
            'y': y,
            't': t,
            'vx': vx,
            'vy': vy,
            'iterations': i}
