'''
Heat advection phenomena - Serial Numpy Vectorized code

Reference : https://www.dropbox.com/s/pkhxlfs1tuftn4w/L5_PhysBased_Unsteady_CHAdvection.pdf
            Slides : 5.39, 5.41, 5.37
'''

import numpy as np
import user_func as uf
from mpi4py import MPI
import sys
import matplotlib.pyplot as plt
import pylab

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

'''

Scheme to approximate advection phenomena
                 1 . First Order Upwind, FOU
                 2 . Second Order Upwind , SOU
                 3 . QUICK
Enter below appropriate serial number of scheme to process :

'''
#scheme = int(raw_input())
scheme = int(sys.argv[1])
assert scheme<4 and scheme>0, 'Enter any one of the three choices in range [1,2,3]'

# Length and height of the problem domain
L, H = 1.0, 1.0
# Maximum number of grid points in L and H
imax, jmax = 162, 162
if rank == 0:
    assert (jmax-2)%size == 0, 'Number of processors doesnt divide the domain evenly'
local_jmax = (jmax-2)/size

# Height and width of each interior Control Volume (CV)
dx = L/(imax-2)
dy = H/(jmax-2)

# Material properties of square domain
rho = 1000.0 # Density $\frac{Kg}{m^3}$
cp  = 4180.0 # Specific heat capacity $\frac{W}{mK}$

# Flow properties
u, v = 1.0, 1.0 # Velocity in x and y direction $\frac{m}{s}

#epsilon = 1e-6; # Convergence criteria
maxiter = 10000

# Implementation of initial and boundary temperature 
t_initial = 50.0
t_left = 100.0
t_top = 100.0 if rank == 0 else t_initial
t_bottom = 0.0 if rank == size-1 else t_initial
t_right = 0.0

t_global = np.zeros((jmax-2,imax)) if rank == 0 else None # temperature profile array
t = np.zeros((local_jmax+2,imax)) + t_initial
t[:,0] = t_left
t[-1,:] = t_bottom 
t[:,-1] = t_right
t[0,:] = t_top
#t_old = t.copy() # Copy of initialized temperature profile array

#Neighbour information
down = rank + 1 
up = rank - 1 

# Computation of time step
dt = (0.2*dy)/np.abs(u)

# Initialization of advection variables
# Advection across CV boundary
t_x = np.zeros((local_jmax,imax-1))
t_y = np.zeros((local_jmax+1,imax-2))
t_x[:,0] = t[1:-1,0]
t_x[:,-1] = t[1:-1,-1]
t_y[0,:] = t[0,1:-1] if rank == 0 else 0.0
t_y[-1,:] = t[-1,1:-1] if rank == size-1 else 0.0
# Net advection flux in the interior CV's
Q = np.zeros((local_jmax,imax-2))

''' In this problem, since velocity is considered to be constant and
uniform through out the domain, mass flow rate doesn't change at all'''
# Calculation of mass flow rate through unit area in x and y directions
mx = rho*u
my = rho*v

# Calculation of width of each CV in x and y directions
x_width = np.zeros((local_jmax,imax))
y_width = np.zeros((local_jmax+2,imax-2))
x_width[:,1:-1] = dx
y_width[:,:] = dy
y_width[0,:] = 0.0 if rank == 0 else dy
y_width[-1,:] = 0.0 if rank == size-1 else dy

'''Calculation of weightage values for interpolation or extrapolation of temperature 
at the interior faces of CV
Note : Here mass flow rate is considered to be positive and uniform through out the domain, So from 
user_func only positive weightages are called to the main function'''
wpx1,wpx2 = uf.weightx(x_width,scheme)
wpy1t,wpy2t = uf.weighty(y_width,scheme)
jtemp,itemp = np.shape(wpy1t)
if rank > 0 and rank < size-1:
    wpy1 = np.zeros((jtemp+2,itemp))
    wpy1[0:-2,:] = wpy1t
    wpy1[-2:,:] = wpy1t[0:2,:]
    wpy2 = np.zeros((jtemp+2,itemp))
    wpy2[0:-2,:] = wpy2t
    wpy2[-2:,:] = wpy2t[0:2,:]
elif rank == 0:
    wpy1 = np.zeros((jtemp+1,itemp))
    wpy1[0:-1,:] = wpy1t
    wpy1[-1,:] = wpy1t[-1,:]
    wpy2 = np.zeros((jtemp+1,itemp))
    wpy2[0:-1,:] = wpy2t
    wpy2[-1,:] = wpy2t[-1,:]
elif rank == size-1:
    wpy1 = np.zeros((jtemp+1,itemp))
    wpy1[1:,:] = wpy1t
    wpy1[0,:] = wpy1t[1,:]
    wpy2 = np.zeros((jtemp+1,itemp))
    wpy2[1:,:] = wpy2t
    wpy2[0,:] = wpy2t[1,:]
    

##if rank == size-5:
##    print wpx1,'\n',wpx2,'\n',wpx3,'\n',wpy1,'\n',wpy2,'\n',wpy3


constant_a = dt/(rho*cp*dx*dy)
iteration = 0

sendbuf_d = np.zeros(imax)
sendbuf_u = np.zeros(imax)
recvbuf_u = np.zeros(imax)
recvbuf_d = np.zeros(imax)

while iteration < maxiter:
    iteration += 1
    #Temperature interpolated or extrapolated in the interior CV faces according to advection scheme
    if rank < size-1:
        sendbuf_d = t[-2,:]
        comm.Send(sendbuf_d, dest = down)
    if rank > 0:
        sendbuf_u = t[1,:]
        comm.Recv(recvbuf_u, source = up)
        t[0,:] = recvbuf_u
        comm.Send(sendbuf_u, dest = up)
    if rank < size-1:
        comm.Recv(recvbuf_d, source = down)
        t[-1,:] = recvbuf_d
    t_x[:,1:-1] = wpx1*t[1:-1,2:-1] + wpx2*t[1:-1,1:-2]
    if rank > 0 and rank < size-1:
        t_y = wpy1*t[0:-1,1:-1] + wpy2*t[1:,1:-1]
    if rank == 0:
        t_y[1:,:] = wpy1*t[1:-1,1:-1] + wpy2*t[2:,1:-1]
    if rank == size-1:
        t_y[0:-1,:] = wpy1*t[0:-2,1:-1] + wpy2*t[1:-1,1:-1]
    adv_x = mx*cp*dy*t_x
    adv_y = my*cp*dx*t_y
    q_adv = (adv_x[:,1:]-adv_x[:,0:-1]) + (adv_y[0:-1,:]-adv_y[1:,:])
    t[1:-1,1:-1] = t[1:-1,1:-1] - constant_a*q_adv
T = t[1:-1,:]
comm.Gather(T,t_global,0)
if rank == 0:
    plt.imshow(t_global)
    pylab.show()
    #print t_global
    
