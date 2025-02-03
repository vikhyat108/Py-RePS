import os
from scipy.integrate import quad
import numpy as np
import math

# Define the range and number of points
y_min = 0.0001
y_max = 1000.0
num_points = 100000  # Adjust the number of points as needed

# Generate the array with logarithmic spacing
y1 = np.logspace(np.log10(y_min), np.log10(y_max), num_points)

# Print the generated array
#print(ytab)

n = len(y1)

def F(x):
    return x*x*np.sqrt(x*x + y*y) / (1 + np.exp(x))

def G(x):
    return x*x/(np.sqrt(x*x + y*y)*(1 + np.exp(x)))

FF = [0.]*n
GG = [0.]*n

file_1 = open("FF_GG.txt","w")

for i in range(0,n):
	#f = F(x,y[i])
	y = y1[i] 
	FF[i], error = quad(F, 0, np.inf)
	#print(f"Result of the integral: {result}")
	
	result, error = quad(G, 0, np.inf)
	GG[i] = y*result
	file_1.write('	'+str(y)+'	'+str(FF[i])+'	'+str(GG[i])+'\n')
	#print(f"Result of the integral: {GG}")




