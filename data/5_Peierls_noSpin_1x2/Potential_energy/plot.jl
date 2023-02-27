using PyPlot

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05

function U(X,w,t1,cm,cp)
	return 0.5* w^2 .* X .^2 - sqrt(2*w)*t1*(cm-cp) .* X
end

Xrng = Array(-30:0.01:30)
w = 0.1
t1 = 0.2
U0 = U(Xrng,w,t1,0.0,0.0)
Up = U(Xrng,w,t1,0.0,1.0)
Um = U(Xrng,w,t1,1.0,0.0)
Upm = U(Xrng,w,t1,0.5,0.5)

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

Us = [U0,Up,Um,Upm]
labels = [L"$\vert c_+\vert^2 = \vert c_-\vert^2 = 0$",
	  L"$\vert c_+\vert^2 = 1,\ \vert c_-\vert^2 = 0$",
	  L"$\vert c_+\vert^2 = 0,\ \vert c_-\vert^2 = 1$",
	  L"$\vert c_+\vert^2 = \vert c_-\vert^2 = \frac{1}{2}$"]

for (u,lab,i) in zip(Us,labels,1:length(Us))
	if i == 4
		lst = "dashed"
	else
		lst = "solid"
	end
	ax.plot(Xrng,u,linewidth=lw,label=lab,linestyle=lst)
end
ax.legend(loc="upper center")
ax.set_xlabel("Vibrational displacement "*L"$X$")
ax.set_ylabel(L"$U(X)\,[\mathrm{eV}]$")
plt.savefig("L2_Peierls_MF_potential.pdf")

