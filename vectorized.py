'''
Heat advection phenomena - Serial Numpy Vectorized code

Reference : https://www.dropbox.com/s/pkhxlfs1tuftn4w/L5_PhysBased_Unsteady_CHAdvection.pdf
            Slides : 5.39, 5.41, 5.37
'''

import numpy as np
import user_func as uf

print '''

Scheme to approximate advection phenomena
                 1 . First Order Upwind, FOU
                 2 . Second Order Upwind , SOU
                 3 . QUICK
Enter below appropriate serial number of scheme to process :

'''
scheme = int(raw_input())
assert scheme<4 and scheme>0, 'Enter any one of the three choices in range [1,2,3]'

# Length and height of the problem domain
L, H = 1.0, 1.0
# Maximum number of grid points in L and H
imax, jmax = 10, 10
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
adv_x = np.zeros((jmax-2,imax-1))
adv_y = np.zeros((jmax-1,imax-2))
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

wpx1,wpx2,wpx3 = uf.weightx(x_width,scheme)
wpy1,wpy2,wpy3 = uf.weighty(y_width,scheme)

constant_a = dt/(rho*cp*dx*dy)
flag = 1

while flag > epsilon:
    #advection flux calculation
    adv_x = mx*(wpx1*t[1:-1,2:-1] + wpx2*t[])


#print t

