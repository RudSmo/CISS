using LinearAlgebra

L = 6;
H = zeros(L,L);

for i in 1:L-1
        H[i,i+1] = -1.0
        H[i+1,i] = -1.0
end

ev = eigvecs(H)

lL = 2
lC = 2

HL = H[1:lL,1:lL]
HC = H[lL+1:lL+lC,lL+1:lL+lC]
HR = H[L-lL+1:L,L-lL+1:L]


el,vl = eigen(HL)
ec,vc = eigen(HC)
er,vr = eigen(HR)

TH = zeros(L,L)
for i in 1:lL
	TH[i,i] = el[i]
end
for i in lL+1:lL+lC
	TH[i,i] = ec[i-lL]
end
for i in L-lL+1:L
	TH[i,i] = er[i-L+lL]
end

#display(TH)

Cl = H[1:lL,1:lL]
Cc = H[lL+1:lL+lC,lL+1:lL+lC]
Cr = H[L-lL+1:L-lL+1]

OAB1 = H[1:lL,lL+1:lL+lC]
OBA1 = H[lL+1:lL+lC,1:lL]
OAB2 = H[lL+1:lL+lC,L-lL+1:L]
OBA2 = H[L-lL+1:L,lL+1:lL+lC]

X1 = vl * OAB1 * vc'
X2 = vc * OBA1 * vl'

Y1 = vr * OAB2 * vc'
Y2 = vc * OBA2 * vr' 


for i in 1:lL
	for j in lL+1:lL+lC
		TH[i,j] = X1[i,j-lL]
		TH[j,i] = X2[j-lL,i]
	end
end
for i in L-lL+1:L
	for j in lL+1:lL+lC
		TH[i,j] = Y1[i-L+lL,j-lL]
		TH[j,i] = Y2[j-lL,i-L+lL]
	end
end
display(TH)
print("\n")
println(eigvals(TH))
println(eigvals(H))

