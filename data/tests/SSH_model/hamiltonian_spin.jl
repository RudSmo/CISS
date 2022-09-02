using LinearAlgebra
using PyPlot


function SSHSpin(L,t1u,t2u,t1d,t2d)
	Hu = zeros(L,L)
	Hd = zeros(L,L)
	for i in 1:2:L-1
		Hu[i,i+1] = -t1u
		Hu[i+1,i] = -t1u
		if i < L-2
			Hu[i+1,i+2] = -t2u
			Hu[i+2,i+1] = -t2u
		end
	end
	
	for i in 1:2:L-1
		Hd[i,i+1] = -t1d
		Hd[i+1,i] = -t1d
		if i < L-2
                        Hd[i+1,i+2] = -t2d
                        Hd[i+2,i+1] = -t2d
                end
	end
	return Hu,Hd
end


L = 70
t1u = 1.0
t2u = 2.0

t1d = 2.0
t2d = 1.0

Hu,Hd = SSHSpin(L,t1u,t2u,t1d,t2d)

ewU,ewD = eigvals(Hu),eigvals(Hd)

#plt.imshow(H,aspect="auto")

#plt.plot(ewU,color="red")
#plt.plot(ewD,color="blue")

#spin up
plt.plot(ewU[1:Int(L/2)],color="red")
plt.plot(Array(0:-1:-Int(L/2)+1),ewU[1:Int(L/2)],color="red")

plt.plot(ewU[end:-1:Int(L/2)+1],color="red")
plt.plot(Array(0:-1:-Int(L/2)+1),ewU[end:-1:Int(L/2)+1],color="red")

#spin down
plt.plot(ewD[1:Int(L/2)],color="blue")
plt.plot(Array(0:-1:-Int(L/2)+1),ewD[1:Int(L/2)],color="blue")

plt.plot(ewD[end:-1:Int(L/2)+1],color="blue")
plt.plot(Array(0:-1:-Int(L/2)+1),ewD[end:-1:Int(L/2)+1],color="blue")

plt.show()

