using LinearAlgebra
using PyPlot


function Ribbon(M,tAA,tAB,tBB)
	L = 6+(M-1)*4
	H = zeros(L,L)*(1+0*1im)
	for i in 1:L-1
		H[i,i+1] = -tAB
		H[i+1,i] = -tAB
	end
	n = 1
	for i in 1:2:Int(L/2)-1
		k = L-1-(n-1)*4 
		H[i,i+k] = -tAB
		H[i+k,i] = -tAB
		n += 1
	end
	for i in 1:2:Int(L/2)-1
		H[i,i+2] = -tAA
		H[i+2,i] = -tAA
	end
	for i in L:-2:Int(L/2)+2
		H[i,i-2] = -tBB
		H[i-2,i] = -tBB
	end
	for i in 2:2:Int(L/2)-1
		j1 = L-(i-2)
		j2 = j1-2
		H[i,j1] = -tBB
		H[j1,i] = -tBB
		H[i,j2] = -tBB
		H[j2,i] = -tBB
	end
	for i in L-1:-2:Int(L/2)+1
		j1 = L-i
		j2 = j1+2
		H[i,j1] = -tAA
		H[j1,i] = -tAA
		H[i,j2] = -tAA
		H[j2,i] = -tAA
	end
	return H
end	

M = 80
L = 6+(M-1)*4


tAB = 1.0

tAAs = Array(0.0:0.1:4.0)

ew = zeros(length(tAAs),L)

for (i,tAA) in enumerate(tAAs)
	tAA = 1im*tAA
	tBB = conj(tAA)
	H = Ribbon(M,tAA,tAB,tBB)
	e = real.(eigvals(H))
	for k in 1:L
		ew[i,k] = e[k]
	end
end

#ew = eigvals(H)
#plt.plot(ew)
#plt.imshow(ew[end:-1:1,:],aspect="auto",extent=[1,L,tAAs[1],tAAs[end]])
plt.plot(ew[3,:])
plt.show()

