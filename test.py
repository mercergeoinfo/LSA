#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plot


def func(x,a,b,c):
    return a*x**2 + b*x + c #Refer [1]

x = np.linspace(0,4,50)
y = func(x,2.6,2,3) + 4*np.random.normal(size=len(x)) #Refer [2]


coeff, var_matrix = curve_fit(func,x,y)
variance = np.diagonal(var_matrix) #Refer [3]

SE = np.sqrt(variance) #Refer [4]

#======Making a dictionary to print results========
results = {'a':[coeff[0],SE[0]],'b':[coeff[1],SE[1]],'c':[coeff[2],SE[2]]}

print "Coeff\tValue\t\tError"
for v,c in results.iteritems():
    print v,"\t",c[0],"\t",c[1]
#========End Results Printing=================

y2 = func(x,coeff[0],coeff[1],coeff[2]) #Saves the y values for the fitted model

plot.plot(x,y)
plot.plot(x,y2)

plot.show()