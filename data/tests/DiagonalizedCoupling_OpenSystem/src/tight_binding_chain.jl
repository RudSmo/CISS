using LinearAlgebra

function swapcol!(x,ii,jj)
	for k in axes(x,1) # edited according to next answer
		for (i,j) in zip(ii,jj)
			idata = x[k,i]
			x[k,i] = x[k,j]
			x[k,j] = idata
		end
	end
end

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

for i in L-l:L
	for j in L-l:L
		Hr[i,j] = H[i,j]
	end
end
HR = H[end-l:end,end-l:end];

Hc = zeros(L,L)
for i in l+1:L-l
	for j in l+1:L-l
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
for i in L-l:L-l+1
	for j in L-l:L-l+1
		Hrc[i,j] = H[i,j]
	end
end
HRC = H[l:l+1,l:l+1]


# Diagonalize L, R

el,vl = eigen(HL)
er,vr = eigen(HR)

# Diagonalize C

ec,vc = eigen(HC)


# Transform LC,RC
U = eigvecs(Hlc)

V = eigvecs(Hc)
swapcol!(V,[1,2,3,4,5,11,12,13,14,15],[6,7,8,9,10,16,17,18,19,20])
THlc = (U * V)' * Hlc * U * V
display(THlc)

G = eigvecs(Hrc)
swapcol!(G,[1,20],[15,19])
swapcol!(G,[19],[18])
swapcol!(G,[18],[17])
swapcol!(G,[17],[16])

THrc = (G * V)' * Hrc * G * V
display(THrc)

TH = V + U + G * THlc + THrc


#display(TH)

