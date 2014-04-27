'''
Heat advection phenomena - Serial Numpy Vectorized code

Reference : https://www.dropbox.com/s/pkhxlfs1tuftn4w/L5_PhysBased_Unsteady_CHAdvection.pdf
            Slides : 5.39, 5.41, 5.37
'''

import numpy as np
from time import time
#import sys
import matplotlib.pyplot as plt
import pylab

T1 = time()
# Length and height of the problem domain
L, H = 1000.0,1000.0
# Maximum number of grid points in L and H
imax, jmax = 302,302
# Height and width of each interior Control Volume (CV)
dx = L/(imax-2)
dy = H/(jmax-2)

# Material properties of square domain
rho = 1000.0 # Density $\frac{Kg}{m^3}$
cp  = 4180.0 # Specific heat capacity $\frac{W}{mK}$

# Flow properties
u, v = 1.0, 1.0 # Velocity in x and y direction $\frac{m}{s}

maxiter = 20000 # Number of iterations to converge

# Implementation of initial and boundary temperature 
t_initial = 50.0
t_left = 100.0
t_top = 100.0
t_bottom = 0.0
t_right = 0.0

t = np.zeros((jmax,imax)) + t_initial # temperature profile array
t[:,0] = t_left
t[-1,:] = t_bottom
t[:,-1] = t_right
t[0,:] = t_top

# Computation of time step
dt = (0.1*dx)/np.abs(u)

# Initialization of advection variables
# Advection across CV boundary
t_x = np.zeros((jmax-2,imax-1))
t_y = np.zeros((jmax-1,imax-2))
t_x[:,0] = t[1:-1,0]
t_x[:,-1] = t[1:-1,-1]
t_y[0,:] = t[0,1:-1]
t_y[-1,:] = t[-1,1:-1]


#In this problem, since velocity is considered to be constant and
#uniform through out the domain, mass flow rate doesn't change at all
# Calculation of mass flow rate through unit area in x and y directions
mx = rho*u
my = rho*v

#Creation of weight matrix for the approximation of temperature at the CV face
#For FOU all weights are having value 1
wpx = np.ones((jmax-2,imax-3))
wpy = np.ones((jmax-3,imax-2))

constant_a = dt/(rho*cp*dx*dy)
iterations = 0

while iterations < maxiter:
    iterations += 1
    #Temperature interpolated or extrapolated in the interior CV faces according to advection scheme
    t_x[:,1:-1] = wpx*t[1:-1,1:-2]
    t_y[1:-1,:] = wpy*t[2:-1,1:-1]
    adv_x = mx*cp*dy*t_x
    adv_y = my*cp*dx*t_y
    q_adv = (adv_x[:,1:]-adv_x[:,0:-1]) + (adv_y[0:-1,:]-adv_y[1:,:])
    t[1:-1,1:-1] = t[1:-1,1:-1] - constant_a*q_adv

T2 = time()
print T2-T1
plt.imshow(t)
pylab.show()

