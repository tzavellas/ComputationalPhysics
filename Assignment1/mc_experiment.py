# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:01:02 2017

@author: Anastasios Tzavellas
"""
import numpy as np
import matplotlib.pyplot as plt

def period(L):
    """
    Calculates the theoretical period of a
    pendulum with length L. It assumes that
    the actual value of g is 9.8m/s2
    
    Args:
        L: The length of the pendulum
    
    Returns:
        The period of the pendulum in sec
    """
    g = 9.8
    return 2*np.pi*np.sqrt(L/g)

def Gravity(L, T, dL, dT):
    """
    Calculates the mean value of g and its error
    given a set of measurements of length, period
    and their respective errors
    
    Args:
        L: The length of the pendulum
        T: The period of the pendulum
        dL: The error in length measurement
        dT: The error in period measurement
    
    Returns:
        A dictionary with 2 arrays: the mean
        value of g and the error that value
    """
    g = np.power(2*np.pi, 2) * L /np.power(T, 2)
    dg = np.power(2*np.pi, 2) * np.sqrt( np.power(dL/(T*T), 2) + np.power(2*L*dT/np.power(T, 3), 2) )
    return {'g':g, 'dg':dg}

def experiment(gTheory, L, deltaL, deltaT, N):
    """
    Performs a g-measurement experiment
    
    Args:
        gTheory: Actual value of g
        L: The lengths of the pendulum at which the experiment takes place
        deltaL: The error in length measurement
        deltaT: The error in period measurement
        N: The number of period measurements
    
    Returns:
        A dictionary with the mean value of g for each length,
        the g-error for each length, the length values and the 
        mean values of the measured period
    """
    g = np.zeros(L.size)                        # Store g measurements
    dg = np.zeros(L.size)                       # Store dg errors
    Tmean = np.zeros(L.size)                    # Store mean value of measured periods
    for i, Li in enumerate(L):
        Lii = np.random.normal(Li, deltaL, N)
        Lmean = np.sum(Lii) / N
        T = period(Lmean)                          # Theoretical Values of T at length Li, for a given g
        Tmeas = np.random.normal(T, deltaT, N)  # Measurements of Period, which is a random variable
                                                # with mean value equal to the theoretical value of
                                                # T and variance deltaT
        Tmean[i] = np.sum(Tmeas) / N            # Mean value of measurements
        dT = np.sqrt( np.sum(Tmeas*Tmeas) / N - np.power(Tmean[i], 2) ) / np.sqrt(N - 1) # Statistical Error
        gStruct = Gravity(Li, Tmean[i], deltaL, dT)   # Value of g with error based on measurements
        g[i] = gStruct['g']
        dg[i] = gStruct['dg']
    return {'g':g, 'dg':dg, 'L':L, 'T':Tmean}

def fit(experiment):
    """
    Performs Least Square Fit on the given experiment
    
    Args:
        experiment: The experiment to perform LSF
    
    Returns:
        A dictionary with the LSF value of g, the
        LSF coefficients, and the values used for the fit
    """
    y = np.power(experiment['T'], 2)
    x = experiment['L']
    A, B = np.polyfit(x, y, 1)              # Least squares fit with x->L and y->T^2
    g = np.power(2*np.pi, 2)/A              # Coefficient A gives g: A = (2*pi)^2 / g
    return {'g':g, 'A':A, 'B':B, 'x':x, 'y':y}


gTheory = 9.8
L = np.array([0.2, 0.4, 0.8, 1.0])
deltaL = 0.01
deltaT = 0.1

print("Experiment with deltaL={:.2f}".format(deltaL))
experiment1 = experiment(gTheory, L, deltaL, deltaT, 20)    # Perform experiment 1 with deltaL = 0.01
print("Mean Value Method")
for i in range(experiment1['g'].size):
    print("L={:.1f} : g = {:.2f} +- {:.2f}".format(L[i], experiment1['g'][i], experiment1['dg'][i]))
print("Least Squares Fit Method")
lsq = fit(experiment1)      # Perform LSF for experiment 1
print("g = {:.2f}".format(lsq['g']))
xn = lsq['x']
yn = np.polyval([lsq['A'], lsq['B']], xn)
plt.plot(xn, yn, 'r', label='$\delta L = 0.01$')        # Plot least square line
plt.plot(lsq['x'], lsq['y'], 'r.')                      # Plot measurements
plt.ylabel('$T^2[sec^2]$')
plt.xlabel('$L[m]$')

deltaL=0.05
print("\nExperiment with systematic error={:.2f}".format(deltaL))
experiment2 = experiment(gTheory, L, deltaL, deltaT, 20)    # Perform experiment 2 with deltaL = 0.05
print("Mean Value Method")
for i in range(experiment2['g'].size):
    print("L={:.1f} : g = {:.2f} +- {:.2f}".format(L[i], experiment2['g'][i], experiment2['dg'][i]))
print("Least Squares Fit Method")
lsq = fit(experiment2)
print("g = {:.2f}".format(lsq['g']))
xn = lsq['x']
yn = np.polyval([lsq['A'], lsq['B']], xn)
plt.plot(xn, yn, 'b', label='$\delta L = 0.05$')        # Plot least square line
plt.plot(lsq['x'], lsq['y'], 'b*')    # Plot measurements
plt.legend()
plt.show()