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

T_l = raw_input("Enter the temperature on the left face : ")

T_t = raw_input("Enter the temperature on the top face : ")

T_b = raw_input("Enter the temperature on the bottom face : ")

for i in range(1,imax-1):
    for j in range(1,jmax-1):

        T(i;j)=T_init; # initial temperature

for i in range(1,imax-1):
    T(i;1)=T_l # west

for i in range(1,imax-1):
    T(i;jmax)=T_l; # north

for j in range(1,jmax-1):
    T(1;j) = T_b; # south

for j in range(1,jmax-1):
    T(imax;j)=T(imax-1;j) # east

##################### calculating mass flow rate through the face ##################

mx = u*rho

my = v*rho

mx_p = max(mx,0);mx_n=min(mx,0);
my_p = max(my,0);my_n=min(my,o);

############# Select the schemes ############################

c = raw_input("Select the schemes to be applied 1 ) FOU 2) SOU 3) QUICK")
