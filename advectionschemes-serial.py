# python code for advection schemes
# serial code
# project me766
# rahul joshi and saran s
# 133106002 / 133106001

##############################################

##### geometrical properties of the plate #############

L = 1.0 # length in meters

H = 1.0 # height in meters

imax = raw_input("Enter the number of grid points in x direction : ");

jmax = raw_input("Enter the number of grid points in y direction : ");

dx = L/(imax-2); # width of the cell in x direction

dy = H/(jmax-2); # width of the cell in y direction

############ fluid properties #######################

rho = 1 # density in kg/m**3

Cp = 0.1 # specific heat in W/mK

k = 0.01 # thermal conductivity

mu = 0.1 # kinematic viscosity

u = 1 # velocity in x direction in m/s

v = 0 # velocity in y direction in m/s

err = 10

eps = 1e-06

################################# time calculation for CFL criterion ######

dt_d = 0.5/((k/(rho*Cp))*((1/dx**2)+(1/dy**2))) # time step stability

dt_c = 1/(max(((abs(u)/dx)+(abs(v)/dy)))) # time step courant

dt = min(dt_d,dt_c);

disp(dt);

Dt = raw_input("Enter the value of the time step: ");

#################################### temperature of the walls and initial BC ###

T = ones([imax,jmax]);

T_init=0 # initial temperature

T_l = raw_input("Enter the temperature")
