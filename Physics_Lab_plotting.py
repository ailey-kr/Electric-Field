#Physics 216
#Plotting code to use the whole year!
#Jaylene Naylor
#September 2015, modified Sept 2017, August 2018
#-------------------------------------------#
#Import packages and libraries needed and give them shortcut names
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def rul3(dA, dB):
    dQ = np.dqrt(dA**2 + dB**2)
    return dQ

#if q = cA^a B^b where c, m, n are constants, then
def rul4(Q, A, dA, a, B, dB, b):
    dQ = Q * np.sqrt(((a*(dA/A))**2) + (b*(dB/B))**2)
    return dQ

def lnbruteforce(A, Aun):
    lnAun = np.abs(np.log(np.average(A)) - np.log(np.average(A) + Aun))
    return lnAun

    
r = np.array([1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7])
V = np.array([3.7, 4.75, 5.35, 5.82, 6.65, 7.1, 7.5, 7.7, 8., 8.25, 8.55, 8.8, 9.13])

r_un = .0005
V_un = .05

#Sigfigs to two places past decimal.
dV = np.diff(V)
dV_un = rule3(V_un, V_un)

dr = np.diff(-1 * r)
dr_un = rule3(r_un, r_un)

Efield = ((-1 * dV)/dr)
Efield_un = (rule4(Efield, 1, -1, dV_un, dr_un, dV, dr)); Efield_unp = np.around(Efield_un, 0)

rmid = (.5 * (r[1:] + r[:-1]))
rmid_un = .5 * (rule3(r_un, r_un))

lnEfield = np.log(Efield); lnEfieldp = np.around(lnEfield, 2)
lnEfield_un = lnbruteforce(Efield, Efield_un); lnEfield_unp = np.around(lnEfield_un, 2)

lnrmid = np.log(rmid); lnrmidp = np.around(lnrmid, 2)
lnrmid_un = lnbruteforce(rmid, rmid_un)
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
print('Change in Voltage:', dV, '+/-', '%.2f' % dV_un, 'volts.\n')
print('E-Field strength', Efield, ' +/-', Efield_unp, 'V/m.\n')
print ('ln Electric Fields:', lnEfieldp, '+/-', lnEfield_unp, '.\n')
print ('ln Average Radii:', lnrmidp, '+/-', '%.2f' % lnrmid_un, '.\n')

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

