# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/ttt/.spyder2/.temp.py
"""

# python code for advection schemes
# serial code
# project me766
# rahul joshi and saran s
# 133106002 / 133106001

##############################################

import numpy as np

import math

import time

##### geometrical properties of the plate #############

L = 1.0 # length in meters

H = 1.0 # height in meters

imax = 12#int(raw_input("Enter the number of grid points in x direction : "));

jmax = 12#int(raw_input("Enter the number of grid points in y direction : "));

dx = L/(imax-2); # width of the cell in x direction

dy = H/(jmax-2); # width of the cell in y direction

############ fluid properties #######################

rho = 1 # density in kg/m**3

Cp = 0.1 # specific heat in W/mK

k = 0.01 # thermal conductivity

mu = 0.1 # kinematic viscosity

u = 1.0 # velocity in x direction in m/s

v = 0.0 # velocity in y direction in m/s

err = 10

eps = 1e-06

################################# time calculation for CFL criterion ######

dt_d = (0.5/((k/(rho*Cp))*((1/dx**2)+(1/dy**2)))) # time step stability

dt_c = (1/(((abs(u)/dx)+(abs(v)/dy)))) # time step courant

dt = min(dt_d,dt_c);

print dt;

Dt = float(raw_input("Enter the value of the time step: "));

alpha=(Dt/(rho*Cp*dx*dy));

#################################### temperature of the walls and initial BC ###

T = np.ones([imax,jmax]);

Tnp = np.ones([imax-1,jmax-1]);

qx=np.zeros([imax,jmax]);

qy=np.zeros([imax,jmax]);

hx=np.zeros([imax,jmax]);

hy=np.zeros([imax,jmax]);

Q_advec=np.zeros([imax,jmax]);

Q_condt=np.zeros([imax,jmax]);

Q_convec=np.zeros([imax,jmax]);

T_init=0 # initial temperature

T_l = int(raw_input("Enter the temperature on the left face : "))

T_t = int(raw_input("Enter the temperature on the top face : "))

T_b = int(raw_input("Enter the temperature on the bottom face : "))

T_r = int(raw_input("Enter the temperature on the right face : "))

for i in range(1,imax-1):
    for j in range(1,jmax-1):

        T[i,j]=T_init; # initial temperature

for i in range(0,imax-1):
    T[i,0]=T_l # west

for i in range(0,imax-1):
    T[i,jmax-1]=T_b; # north

for j in range(0,jmax-1):
    T[0,j] = T_t; # south

for j in range(0,jmax-1):
    T[imax-1,j]=T_r # east

# mass flow in the x and y direction

mx = u*rho

my = v*rho

mx_p=max(mx,0)

mx_n=min(mx,0)

my_p=max(my,0)

my_n=min(my,0)
# select the schemes

c=int(raw_input("Select the scheme to be applied 1) FOU 2) SOU 3) QUICK "))

if c == 1:
    w1,w2,w3=0,1,0
    print (w1,w2,w3);
    
    
elif c==2:
    w1,w2,w3=0,1.5,-0.5
    print (w1,w2,w3);
    
    
else:
    w1,w2,w3=3/8,6/8,-1/8
    print (w1,w2,w3);
    
# carry out the calculations

iter = 0;

while err > eps: 
    
    iter=iter+1;
    
    if iter>1:
        T=Tnp;
    
    #conduction flux in X direction
    for i in range(0,imax):
        for j in range(1,jmax-1):
            qx[i,j]=-(k/dy)*(T[i+1,j]-T[i,j]);
            
            
    #conduction flux in Y direction
    for i in range(1,imax-1):
        for j in range(0,jmax-1):
            qy[i,j]=-(k/dy)*(T[i,j+1]-T[i,j]);
    
    
 #enthalpy in X direction
    
    for j in range(1,jmax-1):
        for i in range(0,imax-1):
            if i==0:
                hx[i,j]=Cp*mx*T[i,j];
                
            elif i==imax-2:
                hx[i,j]=Cp*mx*T[i+1,j];            
            else:
                
                Tn_p=w1*T[i+1,j]+w2*T[i,j]+w3*T[i-1,j];
                Te_n=w1*T[i,j]+w2*T[i+1,j]+w3*T[i+2,j];
                hx[i,j]=Cp*(mx_p*Tn_p + mx_n*Te_n);
            
    
    # heat flux in the Y direction
    
    for j in range(0,jmax-1):
        for i in range(1,imax-1):
            if j==0:
                hy[(i,j)]=Cp*my*T[(i,j)];               
            elif j==jmax-2:
                hy[(i,j)]=Cp*my*T[(i,j+1)];            
            else:
                Tn_p=w1*T[(i,j+1)]+w2*T[(i,j)]+w3*T[(i,j-1)];
                Tn_n=w1*T[(i,j)]+w2*T[(i,j+1)]+w3*T[(i,j+2)]; 
                hy[(i,j)]=Cp*(my_p*Tn_p + my_n*Tn_n);
                
                
    #advection 
    temp=0;
    for i in range(1,imax-1):
        for j in range(1,jmax-1):
            
            Q_advec[(i,j)]=(hx[(i,j)]-hx[(i-1,j)])*dy+(hy[(i,j)]-hy[(i,j-1)])*dx
            
            Q_condt[(i,j)]=(qx[(i-1,j)]-qx[(i,j)])*dy+(qy[(i,j-1)]-qy[(i,j)])*dx
            
            Q_convec[(i,j)]=(Q_condt[(i,j)]-Q_advec[(i,j)]);
            
            Tnp[(i,j)]=T[(i,j)]+alpha*Q_convec[(i,j)];
            
            temp = temp + (Tnp[(i,j)]-T[(i,j)])*(Tnp[(i,j)]-T[(i,j)]);
            
    
    err = abs(math.sqrt(temp/((imax-1)*(jmax-1))));

print ("Iteration = %5d RMS Error (err) = %8.4e\n", iter , err)

 




