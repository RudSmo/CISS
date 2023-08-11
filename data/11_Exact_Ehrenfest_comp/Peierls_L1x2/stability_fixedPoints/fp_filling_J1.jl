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


J(t0,t1,X,c,m) = [0 -1im*(t0-t1*X) -1im*t1*(conj(c)-c) 0;
		   -2*1im*(t0-t1*X) 0 -1im*t1*m 0;
		   0 0 0 1;
		   0 -t1 -1 0]

t0 = 0.04/0.002
t1 = 0.1*t0 # is fixed to 0.004
c1 = 0.1+0.01*1im # imaginary part of coherence on the order of 0.01-0.1
X1 = t0/t1


px = 0.12
py = 0.55
sx = 0.85
sy = 0.42
dy = 0.04

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])

m_rng = -1.0:0.005:1.0 # deviation from half filling, depending on V

ew1 = zeros(length(m_rng),4) * (1+0*1im);
for m_i in 1:length(m_rng)
   ew1[m_i,:] = eigvals(J(t0,t1,X1,c1,m_rng[m_i]))
end


rew = real.(ew1)
ax1.set_ylabel("real part")
ax1.plot((1 .- m_rng) ./ 2,rew,linewidth=lw)
ax1.set_xticks([])


iew = imag.(ew1)
for i in 1:length(iew[:,1])
	iew[i,:] = sort(iew[i,:])
end
ax2.set_ylabel("imaginary part")
ax2.set_xlabel("filling "*L"$n$")
ax2.plot((1 .- m_rng) ./ 2,iew,linewidth=lw)
plt.show()
#plt.savefig("J1_eigenvalues_mVar_t1_$(t1*0.002)_c1_$(abs(c1)).pdf")

