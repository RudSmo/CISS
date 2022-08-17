using LinearAlgebra

L = 20; 
H = zeros(L,L);

for i in 1:L-1
	H[i,i+1] = -1.0
	H[i+1,i] = -1.0
end
l = 5; 

# Split into parts

Hl = zeros(L,L)
for i in 1:l
	for j in 1:l
		Hl[i,j] = H[i,j]
	end
end
HL = H[1:l,1:l];

Hr = zeros(L,L)
for i in end-l:end
	for j in end-l:end
		Hr[i,j] = H[i,j]
	end
end
HR = H[end-l:end,end-l:end];

Hc = zeros(L,L)
for i in l+1:end-l
	for j in l+1:end-l
		Hc[i,j] = H[i,j]
	end
end
HC = H[l+1:end-l,l+1:end-l]

Hlc = zeros(L,L)
for i in l:l+1
	for j in l:l+1
		Hlc[i,j] = H[i,j]
	end
end
HLC = H[l:l+1,l:l+1]

Hrc = zeros(L,L)
for i in l:l+1
	for j in l:l+1
		Hrc[i,j] = H[i,j]
	end
end
HRC = H[l:l+1,l:l+1]


# Diagonalize L, R

el,vl = eigen(HL)
er,vr = eigen(HR)

# Diagonalize C

ec,vc = eigen(HC)
Ec,VC = eigen(Hc)
println(ec)
println(Ec)
# Transform LC,RC

# First, create properly shaped matrices

#uhlcu = 

