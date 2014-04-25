'''
Heat advection phenomena - Serial Numpy Vectorized code

Reference : https://www.dropbox.com/s/pkhxlfs1tuftn4w/L5_PhysBased_Unsteady_CHAdvection.pdf
            Slides : 5.39, 5.41, 5.37
'''

import numpy as np
import user_func as uf
from time import time
import matplotlib.pyplot as plt
import pylab

print '''

Scheme to approximate advection phenomena
                 1 . First Order Upwind, FOU
                 2 . Second Order Upwind , SOU
                 3 . QUICK
Enter below appropriate serial number of scheme to process :

'''
scheme = int(raw_input())
assert scheme<4 and scheme>0, 'Enter any one of the three choices in range [1,2,3]'

t1 = time()
# Length and height of the problem domain
L, H = 1.0, 1.0
# Maximum number of grid points in L and H
imax, jmax = 100, 100
# Height and width of each interior Control Volume (CV)
dx = L/(imax-2)
dy = H/(jmax-2)

# Material properties of square domain
rho = 1000.0 # Density $\frac{Kg}{m^3}$
cp  = 4180.0 # Specific heat capacity $\frac{W}{mK}$

# Flow properties
u, v = 1.0, 1.0 # Velocity in x and y direction $\frac{m}{s}

epsilon = 1e-6; # Convergence criteria

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
t_old = t.copy() # Copy of initialized temperature profile array

# Computation of time step
dt = (0.2*dx)/np.abs(u)

# Initialization of advection variables
# Advection across CV boundary
t_x = np.zeros((jmax-2,imax-1))
t_y = np.zeros((jmax-1,imax-2))
t_x[:,0] = t[1:-1,0]
t_x[:,-1] = t[1:-1,-1]
t_y[0,:] = t[0,1:-1]
t_y[-1,:] = t[-1,1:-1]
# Net advection flux in the interior CV's
Q = np.zeros((jmax-2,imax-2))

''' In this problem, since velocity is considered to be constant and
uniform through out the domain, mass flow rate doesn't change at all'''
# Calculation of mass flow rate through unit area in x and y directions
mx = rho*u
my = rho*v

# Calculation of width of each CV in x and y directions
x_width = np.zeros((jmax-2,imax))
y_width = np.zeros((jmax,imax-2))
x_width[:,1:-1] = dx
y_width[1:-1,:] = dy 

t2 = time()
'''Calculation of weightage values for interpolation or extrapolation of temperature 
at the interior faces of CV
Note : Here mass flow rate is considered to be positive and uniform through out the domain, So from 
user_func only positive weightages are called to the main function'''
wpx1,wpx2 = uf.weightx(x_width,scheme)
wpy1,wpy2 = uf.weighty(y_width,scheme)

t3 = time()

constant_a = dt/(rho*cp*dx*dy)
flag = 1

while flag > epsilon:
    #Temperature interpolated or extrapolated in the interior CV faces according to advection scheme
    t_x[:,1:-1] = wpx1*t[1:-1,2:-1] + wpx2*t[1:-1,1:-2]
    t_y[1:-1,:] = wpy1*t[1:-2,1:-1] + wpy2*t[2:-1,1:-1]
    adv_x = mx*cp*dy*t_x
    adv_y = my*cp*dx*t_y
    q_adv = (adv_x[:,1:]-adv_x[:,0:-1]) + (adv_y[0:-1,:]-adv_y[1:,:])
    t[1:-1,1:-1] = t[1:-1,1:-1] - constant_a*q_adv
    flag = np.sqrt(np.mean((t-t_old)**2))
    t_old = t.copy()

t4 = time()
print t4-t1
plt.imshow(t)
pylab.show()

