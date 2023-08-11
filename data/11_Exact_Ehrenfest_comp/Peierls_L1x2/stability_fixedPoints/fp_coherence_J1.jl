using LinearAlgebra
using PyPlot
plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 1.5
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


J(t0,t1,X,c,m) = [0 -1im*(t0-t1*X) 1im*t1*c 0;
		   -2*1im*(t0-t1*X) 0 -1im*t1*m 0;
		   0 0 0 1;
		   0 -t1 -1 0]

t0 = 0.04/0.002
t1 = 0.1*t0 # is fixed to 0.004
X1 = t0/t1

px = 0.12
py = 0.58
sx = 0.85
sy = 0.39
dy = 0.09

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])

c_rng = -0.2:0.005:0.2
m = 0.0
ew1 = zeros(length(c_rng),4) * (1+0*1im);
for c_i in 1:length(c_rng)
   ew = eigvals(J(t0,t1,X1,c_rng[c_i],m))
   rew = real.(ew)
   iew = imag.(ew)
   ew1[c_i,:] = sort(rew) .+ 1im .* sort(iew)
end


ax1.set_ylabel("real part")
ax1.set_xticks([])
ax1.plot(c_rng,real.(ew1),linewidth=lw)
ax1.set_ylim(-1,1)
ax2.set_ylabel("imaginary part")
ax2.set_xlabel("imaginary part of coherence "*L"$2\mathrm{Im}(c)$")
ax2.plot(c_rng,imag.(ew1),linewidth=lw)
ax2.set_ylim(-1.1,1.1)
plt.savefig("J1_eigenvalues_cVar_t1_$(t1*0.002)_m$(m).pdf")

