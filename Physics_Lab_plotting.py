#Physics 216
#Plotting code to use the whole year!
#Jaylene Naylor
#September 2015, modified Sept 2017, August 2018
#-------------------------------------------#
#Import packages and libraries needed and give them shortcut names
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def rule3(Aun, Bun):
    Qun = np.sqrt(Aun**2 + Bun**2)
    return Qun

def rule4(Q, Ae, Be, Aun, Bun, A, B):
    Qun = np.abs(Q)*np.sqrt(((Ae*Aun/A)**2)+((Be*Bun/B)**2))
    return Qun

def lnbruteforce(A, Aun):
    lnAun = np.abs(np.log(np.average(A)) - np.log(np.average(A) + Aun))
    return lnAun

def rul4(Q, A, dA, a, B, dB, b):
    dQ = Q * np.sqrt(((a*(dA/A))**2) + (b*(dB/B))**2)
    return dQ
    
    
r_array = np.array([1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7])
V_array = np.array([3.7, 4.75, 5.35, 5.82, 6.65, 7.1, 7.5, 7.7, 8., 8.25, 8.55, 8.8, 9.13])

r_un = .0005
V_un = .05

#Sigfigs to two places past decimal.
dV = np.diff(V_array)
dV_un = rule3(V_un, V_un)
print('We had volt changes of', dV, 'with an uncertainty of +/-', '%.2f' % dV_un, 'volts.\n')

#Sigfigs to three places past decimal.
dr = np.diff(-1 * r_array)
dr_un = rule3(r_un, r_un)
print('We had constant changes in radius of', dr[11:], 'with an uncertainty of +/-', '%.3f' % dr_un, 'meters.\n')
#One sigfig.
Efield = ((-1 * dV)/dr)
Efield_un = (rule4(Efield, 1, -1, dV_un, dr_un, dV, dr)); Efield_unp = np.around(Efield_un, 0)
print('The strength of our electric fields (starting at .0675m) were:', Efield, '\n with uncertainties of +/-', Efield_unp, 'volts per meter.\n')

#Three sigfigs?
rmid = (.5 * (r_array[1:] + r_array[:-1]))
rmid_un = .5 * (rule3(r_un, r_un))
print('The averages of our radii were:', rmid, 'with an uncertainty of +/-', '%.4f' % rmid_un, 'meters.\n')

#One sigfig.
lnEfield = np.log(Efield); lnEfieldp = np.around(lnEfield, 2)
lnEfield_un = lnbruteforce(Efield, Efield_un); lnEfield_unp = np.around(lnEfield_un, 2)
print ('ln Electric Fields:', lnEfieldp, '\n with uncertainties of +/-', lnEfield_unp, '.\n')

#Three sigfigs?
lnrmid = np.log(rmid); lnrmidp = np.around(lnrmid, 2)
lnrmid_un = lnbruteforce(rmid, rmid_un)
print ('ln Average Radii:', lnrmidp, 'with an uncertainty of +/-', '%.2f' % lnrmid_un, '.\n')
x = lnrmid   #this should be the array you want to plot on the x axis
y = lnEfield
dy = lnEfield_un  #this should be your error in y array

#----------------------------------------------#
#Don't need to change anything in this section!
 
#Find the intercept and slope, b and m, from Python's polynomial fitting function
b,m=np.polynomial.polynomial.polyfit(x,y,1,w=dy)

#Write the equation for the best fit line based on the slope and intercept
fit = b+m*x

#Calculate the error in slope and intercept
#def Delta(x, dy) is a function, and we will learn how to write our own at a later date. They are very useful!
def Delta(x, dy):
    D = (sum(1/dy**2))*(sum(x**2/dy**2))-(sum(x/dy**2))**2
    return D
 
D=Delta(x, dy)
 
dm = np.sqrt(1/D*sum(1/dy**2)) #error in slope
db = np.sqrt(1/D*sum(x**2/dy**2)) #error in intercept

#Calculate the "goodness of fit" from the linear least squares fitting document
def LLSFD2(x,y,dy):
    N = sum(((y-b-m*x)/dy)**2)
    return N
                      
N = LLSFD2(x,y,dy)

#-----------------------------------------------------------------------#
#Plot data on graph. Plot error bars and place values for slope, error in slope and goodness of fit on the plot using "annotate"
plt.figure(figsize=(15,10))
 
plt.plot(x, fit, color='green', linestyle='--')
plt.scatter(x, y, color='blue', marker='o')
 
#create labels  YOU NEED TO CHANGE THESE!!!
plt.xlabel('ln Radius')
plt.ylabel('ln Electric Field')
plt.title('Electric Field w/ regard to Radius')
 
plt.errorbar(x, y, yerr=dy, xerr=None, fmt="none") #don't need to plot x error bars
 
plt.annotate('Slope (unitless) = {value:.{digits}E}'.format(value=m, digits=2),
             (0.05, 0.9), xycoords='axes fraction')
 
plt.annotate('Error in Slope (ditto) = {value:.{digits}E}'.format(value=dm, digits=1),
             (0.05, 0.85), xycoords='axes fraction')
 
plt.annotate('Goodness of fit = {value:.{digits}E}'.format(value=N, digits=2),
             (0.05, 0.80), xycoords='axes fraction')

plt.show()

