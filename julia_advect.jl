T1 = time()
H,L = 1000.0,1000.0

jmax,imax = 302,302

dx = L/(imax-2)
dy = H/(jmax-2)

rho,cp = 1000.0,4180.0

u,v = 1.0,1.0

maxiter = 20000

t_initial = 50.0
t_left = 100.0
t_top = 100.0
t_bottom = 0.0
t_right = 0.0

t = zeros(jmax,imax)
t[:,1] = t_left
t[end,:] = t_bottom
t[:,end] = t_right
t[1,:] = t_top

dt = (0.1*dx)/u

t_x = zeros(jmax-2,imax-1)
t_y = zeros(jmax-1,imax-2)
t_x[:,1] = t[2:end-1,1]
t_x[:,end] = t[2:end-1,end]
t_y[1,:] = t[1,2:end-1]
t_y[end,:] = t[end,2:end-1]

mx = rho*u
my = rho*v

wpx = ones(jmax-2,imax-3)
wpy = ones(jmax-3,imax-2)

constant_a = dt/(rho*cp*dx*dy)
iterations = 0

while iterations < maxiter
    iterations += 1
    t_x[:,2:end-1] = wpx.*t[2:end-1,2:end-2]
    t_y[2:end-1,:] = wpy.*t[3:end-1,2:end-1]
    adv_x = mx*cp*dy*t_x
    adv_y = my*cp*dx*t_y
    q_adv = (adv_x[:,2:end]-adv_x[:,1:end-1]) + (adv_y[1:end-1,:]-adv_y[2:end,:])
    t[2:end-1,2:end-1] = t[2:end-1,2:end-1] - constant_a*q_adv
end
T2 = time()
println(T2-T1)
