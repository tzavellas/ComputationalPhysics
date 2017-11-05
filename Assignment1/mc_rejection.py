# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 18:09:06 2017

@author: Anastasios Tzavellas
"""

import numpy as np
import matplotlib.pyplot as plt

# ------------------ Start of the program ------------------

N = 10000     # Number of samples
xInterval = np.array([0, 2]) 	# Interval of x values
yInterval = np.array([0, 2]) 	# Interval of y values
area = xInterval[1] * yInterval[1] 	# Total area of Rejection method

x = np.random.uniform(0, xInterval[1], N) 	# Generate N random values in U(0, xInterval[1])
y = np.random.uniform(0, yInterval[1], N) 	# Generate N random values in U(0, yInterval[1])

inside = np.where( (x+y<2) & (x*x+y*y>1) ) 	# inside contains the indices of the x and y samples
													# that lie within  the boundaries of the object
density = 1+ np.power(x[inside], 2) + np.power(y[inside], 2) # density contains the density values at the above indices

sumMass = np.sum(density)						# The sumMass is the sum of the density values
sumSigmaMass = np.sum(density * density) 		# The sumSigmaMass is the sum of the square of the density values

sumXcm = np.sum(x[inside]*density) 	 	 	# The sumXcm is the sum of the x*density values
sumSigmaXcm = np.sum((x[inside]*density) * (x[inside]*density)) # The sumSigmaXcm is the sum of the square of the x*density values

sumYcm = np.sum(y[inside]*density) 			# The sumYcm is the sum of the y*density values
sumSigmaYcm = np.sum((y[inside]*density) * (y[inside]*density)) # The sumSigmaYcm is the sum of the square of the y*density values

mass = area * sumMass / N 						# The mass of the object
errorMass = mass * np.sqrt(1-inside[0].size/N) / np.sqrt(inside[0].size) # The error of the calculation

x_cm = (area * sumXcm/N) / mass 				# The x coordinate of the center of mass of the object
errorXcm = x_cm * np.sqrt(1-inside[0].size/N) / np.sqrt(inside[0].size) # The error of the calculation

y_cm = (area * sumYcm/N) / mass 				# The y coordinate of the center of mass of the object
errorYcm = y_cm * np.sqrt(1-inside[0].size/N) / np.sqrt(inside[0].size) # The error of the calculation

print("Monte Carlo Rejection Method (", N, " samples)\n")
print('mass = {:.3f} +- {:.3f}'.format(mass, errorMass))
print('x_cm = {:.3f} +- {:.3f}'.format(x_cm, errorXcm))
print('y_cm = {:.3f} +- {:.3f}'.format(y_cm, errorYcm))

print("\nAnalytic Evaluation for Comparison\n")
print('mass = 3.48857')
print('x_cm = 0.840841')
print('y_cm = 0.840841')

plt.plot(x[inside], y[inside], 'k.', markersize=2)
plt.show()