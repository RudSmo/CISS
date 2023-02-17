module ClassicalSpins

using LinearAlgebra
export SpinEOM

function normMag(mag,nfac)
	normed = []
	for mg in mag
		nm = norm(mg)
		if nm != 0.0
			push!(normed,nfac*mg/nm)
		else
			push!(normed,[0.,0.,0.])
		end
	end
	return normed
end

function SpinEOM(Jsd,Scd,M,pos,Hext)
	L = []
	for i in 1:length(pos)
		m = Jsd*cross(Scd[pos[i]],M[i])-cross(Hext[i],M[i])
		push!(L,m)
	end
	return L
end

end

