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

def gError(L, T, dL, dT):
    """
    Calculates the error of g 
    estimation given the length, the 
    period and their respective errors
    using the error propagation rule
    
    Args:
        L: A vector of length values
        T: A vector of period values
        dL: The error in length measurement
        dT: The error in period measurement
    
    Returns:
        A vector with g error values
    """
    dg = np.power(2*np.pi, 2) * np.sqrt( np.power(dL/(T*T), 2) + np.power(2*L*dT/np.power(T, 3), 2) )
    return dg

def experiment(gTheory, L, dL, dT, dLsystm = 0):
    """
    Performs a g-measurement experiment
    
    Args:
        gTheory: Actual value of g
        L: A vector of length values of the pendulum
        dL: The error in length measurement
        dT: The error in period measurement
        dLsystm: Systematic error of length measurement, default value 0
    
    Returns:
        A dictionary with the mean values of g, 
        the g-error values and the measured period 
        values, for each length
    """
    L = L + dLsystm                                             # Add systematic error, if it exists
    T = period(L) + np.random.standard_normal(L.size) * dT      # Period measurements with dt=0.1s accuracy
    g = np.power(2*np.pi, 2) * L / np.power(T, 2)               # Indirect g measurement from length and period
    gMean = np.sum(g)/g.size                                    # Mean value of g measurements
    dg = gError(L, T, dL, dT)                                   # g measurement error
    return {'g':g, 'dg':dg, 'L':L, 'T':T, 'gMean':gMean}

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
    A, B = np.polyfit(x, y, 1)              # Least squares fit with x->L and y->T^2: y = Ax + B
    g = np.power(2*np.pi, 2)/A              # Coefficient A gives g: A = (2*pi)^2 / g
    return {'g':g, 'A':A, 'B':B, 'x':x, 'y':y}


gTheory = 9.8
L = np.array([0.2, 0.4, 0.8, 1.0])
dL = 0.01
dT = 0.1

print("Experiment without systematic error")
experiment1 = experiment(gTheory, L, dL, dT)    # Perform experiment 1
print("Mean Value Method")
print("-----------------")
for i in range(experiment1['g'].size):
    print("L={:.1f} : g = {:.2f} +- {:.2f}".format(L[i], experiment1['g'][i], experiment1['dg'][i]))
print("g Mean value: {:.2f}".format(experiment1['gMean']))
print("\nLeast Squares Fit Method")
print("------------------------")
lsq1 = fit(experiment1)                          # Perform LSF for experiment 1
print("g Least Squares: {:.2f}".format(lsq1['g']))
xn = lsq1['x']
yn = np.polyval([lsq1['A'], lsq1['B']], xn)
plt.plot(xn, yn, 'r', label='$\delta L_{system} = 0$')# Plot least square line
plt.errorbar(lsq1['x'], lsq1['y'], xerr=dL, yerr=dT, fmt='r.')              # Plot measurements
plt.ylabel('$T^2[sec^2]$')
plt.xlabel('$L[m]$')

dLsystm = 0.05
print("\n\nExperiment with systematic error={:.2f}".format(dLsystm))
experiment2 = experiment(gTheory, L, dL, dT, dLsystm)   # Perform experiment 2 with dL = 0.05
print("Mean Value Method")
print("-----------------")
for i in range(experiment2['g'].size):
    print("L={:.1f} : g = {:.2f} +- {:.2f}".format(L[i], experiment2['g'][i], experiment2['dg'][i]))
print("g Mean value: {:.2f}".format(experiment2['gMean']))
print("\nLeast Squares Fit Method")
print("-------------------------")
lsq2 = fit(experiment2)
print("g Least Squares: {:.2f}".format(lsq2['g']))
xn = lsq2['x']
yn = np.polyval([lsq2['A'], lsq2['B']], xn)
plt.plot(xn, yn, 'b', label='$\delta L_{system} = 0.05$')        # Plot least square line
plt.errorbar(lsq2['x'], lsq2['y'], xerr=dL, yerr=dT, fmt='b*')                      # Plot measurements
plt.legend()
plt.show()