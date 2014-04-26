'''
Heat advection phenomena - MPI Vectorized code

Reference : https://www.dropbox.com/s/pkhxlfs1tuftn4w/L5_PhysBased_Unsteady_CHAdvection.pdf
            Slides : 5.39, 5.41, 5.37
'''

import numpy as np
from mpi4py import MPI
import sys
import matplotlib.pyplot as plt
import pylab
import time

t1=time.time()

# Setting up the computational world for the program
# Here default communication world is used
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Length and height of the problem domain
L, H = np.float(sys.argv[1]),np.float(sys.argv[1])
# Maximum number of grid points in L and H including boundary CV
imax, jmax = np.int(sys.argv[2]), np.int(sys.argv[2])
if rank == 0:
    assert (jmax-2)%size == 0, 'Number of processors doesnt divide the domain evenly'
# Total domain is sub-divided into smaller blocks along the y or H direction
# Each of these sub-divisions are individually calculated in each processor
# So now each processor has a problem domain of matrix size $local_jmax+2 \times imax$
# which includes 2 extra rows padded to each sub-division to track the neighbours from other processors
local_jmax = (jmax-2)/size

# Height and width of each interior Control Volume (CV)
dx = L/(imax-2)
dy = H/(jmax-2)

# Material properties of square domain
rho = 1000.0 # Density $\frac{Kg}{m^3}$
cp  = 4180.0 # Specific heat capacity $\frac{W}{mK}$

# Flow properties
u, v = 1.0, 1.0 # Velocity in x and y direction $\frac{m}{s}

# Number of iterations to run to get the steady state solution
maxiter = np.int(sys.argv[3])

# Implementation of initial and boundary temperature 
# Top and Bottom boundary is only updated in two process
# Rest all processes are initialised with bottom and top as initial temperature
t_initial = 50.0
t_left = 100.0
t_top = 100.0 if rank == 0 else t_initial
t_bottom = 0.0 if rank == size-1 else t_initial
t_right = 0.0

# t_global gathers temperature information from all working processes for post procesing of data
t_global = np.zeros((jmax-2,imax)) if rank == 0 else None # temperature profile array
# initialising local temperature array for calculations
t = np.zeros((local_jmax+2,imax)) + t_initial
t[:,0] = t_left
t[-1,:] = t_bottom 
t[:,-1] = t_right
t[0,:] = t_top

#Neighbour process information to ensure the continuity for the breakdown process
down = rank + 1 
up = rank - 1 

# Computation of time step
dt = (0.1*dy)/np.abs(u)

# Initialization of advection variables
# Advection across CV boundary in both X and Y direction
t_x = np.zeros((local_jmax,imax-1))
t_y = np.zeros((local_jmax+1,imax-2))
# Boundary temperatures need not be approximated using advection scheme
# Boundary temperature values are directly copied
t_x[:,0] = t[1:-1,0]
t_x[:,-1] = t[1:-1,-1]
t_y[0,:] = t[0,1:-1] if rank == 0 else 0.0
t_y[-1,:] = t[-1,1:-1] if rank == size-1 else 0.0

#In this problem, since velocity is considered to be constant and
#uniform through out the domain, mass flow rate doesn't change at all
# Calculation of mass flow rate through unit area in x and y directions
mx = rho*u
my = rho*v

#Creation of weight matrix for the approximation of temperature at the CV face
#For FOU all weights are having value 1
#For inner sub-divisions of domain the weight matrix size is different from that of boundary sub-divisions 
wpx = np.ones((local_jmax,imax-3))
wpy = np.ones((local_jmax+1,imax-2)) if rank > 0 and rank < size-1 else np.ones((local_jmax,imax-2))
    
constant_a = dt/(rho*cp*dx*dy)
iteration = 0

# Initilising the buffers need to be send and recieved 
# Which ensures the continuity of the problem
sendbuf_d = np.zeros(imax)
sendbuf_u = np.zeros(imax)
recvbuf_u = np.zeros(imax)
recvbuf_d = np.zeros(imax)

while iteration < maxiter:
    iteration += 1

    # Transfer of data between process to ensure the continuity
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

    # Calculation of temperature at CV faces using FOU advection scheme
    t_x[:,1:-1] = wpx*t[1:-1,1:-2]
    if rank > 0 and rank < size-1:
        t_y = wpy*t[1:,1:-1]
    if rank == 0:
        t_y[1:,:] = wpy*t[2:,1:-1]
    if rank == size-1:
        t_y[0:-1,:] = wpy*t[1:-1,1:-1]
    
    # Physical method of calculation of temperature as discussed in reference at the start of code
    adv_x = mx*cp*dy*t_x
    adv_y = my*cp*dx*t_y
    q_adv = (adv_x[:,1:]-adv_x[:,0:-1]) + (adv_y[0:-1,:]-adv_y[1:,:])
    t[1:-1,1:-1] = t[1:-1,1:-1] - constant_a*q_adv

# Collecting the updated result to a local variable which will be send to the root process later
T = t[1:-1,:]
comm.Gather(T,t_global,0)
if rank == 0:
    plt.imshow(t_global)
    pylab.show()

t2=time.time()

t3=t2-t1

print t3
