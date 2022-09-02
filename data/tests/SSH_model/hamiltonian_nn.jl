using LinearAlgebra
using PyPlot

function NNTB(L,t1,t2)
	H = zeros(L,L)
	for i in 1:L-1
		H[i,i+1] = -t1
		H[i+1,i] = -t1
		if i < L-2
			H[i,i+2] = -t2
			H[i+2,i] = -t2
		end
	end
	return H
end


L = 150
t1 = 1.0
t2 = 2.0
H = NNTB(L,t1,t2)
ew = eigvals(H)
#plt.imshow(H,aspect="auto")
plt.plot(ew)

#plt.plot(ew[1:Int(L/2)],color="k")
#plt.plot(Array(0:-1:-Int(length(ew)/2)+1),ew[1:Int(L/2)],color="k")
#
#plt.plot(ew[end:-1:Int(L/2)+1],color="k")
#plt.plot(Array(0:-1:-Int(length(ew)/2)+1),ew[end:-1:Int(L/2)+1],color="k")
plt.show()

