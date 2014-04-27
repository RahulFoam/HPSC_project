T1 = time()
H,L = 1000.0,1000.0

jmax,imax = 12,12

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
for i in 1:jmax
	t[i,1] = t_left
	t[end,i] = t_bottom
	t[i,end] = t_right
	t[1,i] = t_top
end
dt = (0.1*dx)/u

t_x = zeros(jmax-2,imax-1)
t_y = zeros(jmax-1,imax-2)
for i in 1:jmax-2
	t_x[i,1] = t[i+1,1]	
	t_x[i,end] = t[i+1,end]
	t_y[1,i] = t[1,i+1]
	t_y[end,i] = t[end,i+1]
end

mx = rho*u
my = rho*v

wpx = ones(jmax-2,imax-3)
wpy = ones(jmax-3,imax-2)

constant_a = dt/(rho*cp*dx*dy)
q_adv = zeros(jmax-2,imax-2)

for iter in 1:maxiter
        for j in 1:jmax-2
		for i in 1:imax-3
      			t_x[j,i+1] = wpx[j,i]*t[j+1,i+1]
		end
	end
        for j in 1:jmax-3
		for i in 1:imax-2
      			t_y[j+1,i] = wpy[j,i]*t[j+2,i+1]
		end
	end 
    	adv_x = mx*cp*dy*t_x
    	adv_y = my*cp*dx*t_y
	for j in 1:jmax-2
		for i in 1:imax-2
    			q_adv[j,i] = (adv_x[j,i+1]-adv_x[j,i]) + (adv_y[j,i]-adv_y[j+1,i])
    			t[j+1,i+1] = t[j+1,i+1] - constant_a*q_adv[j,i]
		end
	end
end
T2 = time()
println(T2-T1)
println(t)
