include("../../../src/exact_TDSE/hamiltonian.jl")
using .Hamiltonian
using PyPlot
using LinearAlgebra
using SparseArrays
using Combinatorics 
using StatsBase

function Idkron(n::Int64,b::Int64)::Any
    """
    Tensor product of n identity operators 1(b,b).
    Input:
        n ::Int64, number of Identity operators
        b ::Int64, dimension of identity operator
    """
    id = Matrix(I,b,b)
    if n == 0
        return 1.0
    elseif n == 1
        return id
    else
        return kron(id,Idkron(n-1,b))
    end
end
function bop(M::Int64)::Array
    """
    Bosonic annihilation operator with M levels.
    Input:
        M ::Int64, number of levels
    """
    b = zeros(M,M)
    for i in 1:M-1
                b[i,i+1] = sqrt(i)
        end
    return b
end
function nFermi(op,i,num)
	id = [1 0; 0 1.]
	sz = [1 0; 0 -1.]
	A = []
	for j in 1:num
		if j < i
			push!(A,sz)
		elseif j == i
			push!(A,op)
		else
			push!(A,id)
		end
	end
	return A		
end

function nBose(op,i,num,vib)
	id = Matrix(I,vib,vib)
	A = []
	for j in 1:num
		if j == i
			push!(A,op)
		else
			push!(A,id)
		end
	end
	return A
end

function OperatorProduct(A,B)
	C = []
	nA = length(A)
	for i in 1:nA
		push!(C,A[i]*B[i])
	end
	return nkron(C)
end
function VDot(a)
    sx = [0 1; 1 0.]
    sy = [0 -1im; 1im 0.]
    sz = [1 0; 0 -1.]
    return a[1]*sx + a[2]*sy + a[3]*sz
end

function absv(a)
    return sqrt(dot(a,a))
end

function UnitVector(a,c,j,num_laps,num_elec)
    phij = 2*pi*(j-1)*num_laps/(num_elec-1)
    cj = c*(j-1)/(num_elec-1)
    return [a*cos(phij),a*sin(phij),cj]
end
function VSVec(num_elec,a_helix,c_helix,site,num_laps)
    rj = UnitVector(a_helix,c_helix,site,num_laps,num_elec)
    rj1 = UnitVector(a_helix,c_helix,site+1,num_laps,num_elec)
    rj2 = UnitVector(a_helix,c_helix,site+2,num_laps,num_elec)
    rjrj1 = rj-rj1
    rjrj2 = rj-rj2

    dj1 = rjrj1/absv(rjrj1)
    dj2 = rjrj2/absv(rjrj2)
    #println(dj1)
    #println(dj2)
    vj = cross(dj1,dj2)
    return VDot(vj)
end

function Fransson(num_sites::Int64,num_laps::Int64,vibStates::Int64,eps::Array,t0::Float64,tPh::Float64,
                  λ0::Float64,λPh::Float64,wPh::Float64,r_helix::Float64,c_helix::Float64)::Array
    """
    Fransson's Hamiltonian, that is a fermionic Hubbard chain with L electronic sites. Each site is coupled to a harmonic oscillator
    with vibState levels. In the electronic chain, electrons can hop between neighboring sites, a hopping induces vibrational transitions and
    there is a SOC between next-nearest-neighboring sites, aus well as a spin-dependent electron-phonon contribution.
    The electron-phonon (g) is for now a scalar, but can be made energy-dependent at some later
    time. Also, it is not possible to specify the number of bosonic levels lattice-site dependently. Each mode has the same number of
    levels. The coupling is nearest-neighbor sites only. I'm generating the Hamiltonian literally, meaning all tensor products are performed explicitly. This gives the Hamiltonian in
    the full Hilbert space. Certainly, both the construction of the Hamiltonian and not restricting oneself to a subspace is not the
    most optimal way to do things, but it works.

    Input:
        L        :: Int64, lattice sites
        num_laps :: Int64, number of laps in chiral molecule
        M        :: Int64, bosonic levels
        eps      :: Array, on-site energies
        t0       :: Float64, nearest-neighbor hopping
        t        :: Float64, nearest-neighbor hopping with electron-phonon interaction
        λ0       :: Float64, SOC between next-nearest neighbors
        λ        :: Float64, SOC between next-nearest neighbors with electron-phonon interaction
        r_hel    :: Float64, radius of helix
        c_hel    :: Float64, size of helix
    Output:
        H ::Array, the Heisenberg-Holstein Hamiltonian.
    """
    d = 2
    eu,ed = eps
    num = num_sites * num_laps
    dim = d^(2*num)*vibStates^(num-1)
    IdS = Matrix(I,d^(2*num),d^(2*num))
    IdB = Matrix(I,vibStates^(num-1),vibStates^(num-1))
    id = Matrix(I,2,2)
    HS = zeros(dim,dim)*1im
    HB = zeros(dim,dim)
    println(dim)
    HSB = zeros(dim,dim)
    c = [0 1; 0 0.]
    cd = c'
    nd = cd*c

    # Harmonic oscillators
    b = bop(vibStates)
    bd = b'
    for i in 1:num
	bdb = nkron(nBose(bd*b,i,num-1,vibStates))
        HB += wPh*kron(IdS,bdb)
    end
    # On-site
    for i in 1:num
	cdu_i = nFermi(cd,2*i-1,2*num)
	cdd_i = nFermi(cd,2*i,2*num)
	cu_i = nFermi(c,2*i-1,2*num)
	cd_i = nFermi(c,2*i,2*num)

        hu = eu * OperatorProduct(cdu_i,cu_i)
        hd = ed * OperatorProduct(cdd_i,cd_i)
        HS += kron(hu + hd,IdB)
    end
    # Hopping
    for i in 1:num-1
	cdu_i = nFermi(cd,2*i-1,2*num)
	cdd_i = nFermi(cd,2*i,2*num)
	cu_i = nFermi(c,2*(i+1)-1,2*num)
	cd_i = nFermi(c,2*(i+1),2*num)
        hu = -t0* OperatorProduct(cdu_i,cu_i)
        hd = -t0* OperatorProduct(cdd_i,cd_i)
	h_hop = hu+hd + hu'+hd'
        HS += kron(hu+hu'+hd+hd',IdB)
	qv   = nkron(nBose(bd+b,i,num-1,vibStates))
	HSB -= tPh*kron(h_hop,qv)	
    end
    # SOC
    for i in 1:num-2
	cdu_i = nFermi(cd,2*i-1,2*num)
        cdd_i = nFermi(cd,2*i,2*num)
        cu_i = nFermi(c,2*(i+2)-1,2*num)
        cd_i = nFermi(c,2*(i+2),2*num)	
        vs  = VSVec(num_sites,r_helix,c_helix,i,num_laps)
        huu = 1im*λ0*vs[1,1]*OperatorProduct(cdu_i,cu_i)
        hud = 1im*λ0*vs[1,2]*OperatorProduct(cdu_i,cd_i)
        hdu = 1im*λ0*vs[2,1]*OperatorProduct(cdd_i,cu_i)
        hdd = 1im*λ0*vs[2,2]*OperatorProduct(cdd_i,cd_i)
	hsoc = huu + hud + hdu + hdd + huu' + hud' + hdu' + hdd'
        HS += kron(hsoc,IdB)
	qv  = nkron(nBose(bd+b,i,num-1,vibStates))
	HSB += kron(hsoc,qv)
    end
    H = HS+HB+HSB
    return H
end

function nkron(A)
	if length(A) == 1
		return A[1]
	else
		kron(A[1],nkron(A[2:end]))
	end
end
function ElectronBasis(sites)
  	A = [("0,"^Int(sites/2-1)*"0",nkron([[1,0] for i in 1:sites]))]
	for i in 1:sites-1
      		on = ones(i)
      		ze = zeros(sites-i)
      		b = collect(Iterators.flatten([ze,on]))
	      	for bi in collect(myperms(b))[end:-1:1]
	  		str = ""
			for j in 1:2:sites
				k1,k2 = bi[j:j+1]
				if k1 == 0 && k2 == 0
					str *= "0"
				elseif k1 == 1 && k2 == 0
					str *= "↑"
				elseif k1 == 0 && k2 == 1
					str *= "↓"
				elseif k1 == 1 && k2 == 1
					str *= "⇅"
				end
				if j != sites-1
					str *= ","
				end
			end
			op = 1 
			for i in bi
				if i == 0
					op = kron(op,[1,0])
				else
					op = kron(op,[0,1])
				end		
			end
	  		push!(A, (str,op))
      		end
  	end
  	push!(A,("⇅,"^Int(sites/2-1)*"⇅",nkron([[0,1] for i in 1:sites])))
  	return A
end
function VibBasis(N)
	on = ones(1)
	ze = zeros(N-1)
	conf = collect(Iterators.flatten([on,ze]))
	perm_conf = collect(myperms(conf))
	conf_all = []
	for ci in perm_conf
		str = string(findall(x->x==1,ci)[1]-1)
		push!(conf_all,(str,ci))
	end
		
	return conf_all
end
function SingleListKron(A,B,iter,N)
	nA = length(A)
	nB = length(B)
	conf = []
	if iter == 1
		s1 = "|"
		s2 = "⟩|"
		s3 = ""
	elseif iter == N-1
		s1 = ""
                s2 = ","
                s3 = "⟩"
	else
		s1 = ""
		s2 = ","
		s3 = ""
	end
	for a in A
		for b in B
			push!(conf,(s1*a[1]*s2*b[1]*s3,kron(a[2],b[2])))
		end
	end
	return conf
end
function MultiListKron(A)
	nA = length(A)
	B = A[1]
	for i in 2:nA
		B = SingleListKron(B,A[i],i-1,nA)
	end
	return B
end
function trans(x, v::Dict{T, Int}, l) where T
   z = collect(1:l)
   idxs = Vector{Int}[]
   for k in x
       push!(idxs, z[k])
       deleteat!(z, k)
   end
   res = Vector{T}(undef, l)
   for (j, k) in enumerate(keys(v))
       for i in idxs[j]
	   res[i] = k
       end
   end
   res
end
function myperms(x)
   v = countmap(x)
   s = Int[length(x)]
   for (k,y) in v
       l = s[end]-y
       l > 0 && push!(s, l)
   end
   iter = Iterators.product((combinations(1:s[i], vv) for (i, vv) in enumerate(values(v)))...)
   (trans(z, v, length(x)) for z in iter)
end

function GenerateFullBasis(nSites,vibStates)
	state_confs = []
	el_basis = ElectronBasis(2*nSites)	
	push!(state_confs,el_basis)
	for i in 1:nSites-1
		push!(state_confs,VibBasis(vibStates))
	end
	return MultiListKron(state_confs)
end


num_sites = 3
num_laps = 1
N = num_sites * num_laps
vibStates = 10
T0 = 0.04
t0 = T0
E =  -6*t0
eps = [E,E]
tPh = 0.1*t0
λ0 = t0/40
λPh = 0.1 * λ0
wPh = 0.002
r_helix = 1.0
c_helix = 1.0


H = Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,λ0,λPh,wPh,r_helix,c_helix)
fil = open("Hamiltonian_Fransson_$(num_laps)x$(num_sites)_eps_$(E)_w0_$(wPh)_t0_$(T0)_t1_$(tPh)_tso_$(λ0)_tsoPh_$(λPh)_nv$(vibStates).txt","w")
sH = sparse(H)
println(sH)
write(fil,H)


#state_confs = GenerateFullBasis(N,vibStates)
#
#
#ev = eigvecs(H)
#println(size(ev))
#GS = ev[:,1]
#
#overlaps = []
#confstr = []
#for conf in state_confs
#	cstr,cstate = conf
#	ol = cstate' * GS
#	if abs(ol) > 1e-3
#		push!(overlaps,real(ol))
#		push!(confstr,cstr)
#	end
#end
#
#so = sortperm(overlaps)
#overlaps = overlaps[so][end:-1:1]
#confstr = confstr[so][end:-1:1]
#
#for (cstr,ol) in zip(confstr,overlaps)
#	println(cstr," ",ol)
#end
#
#
#
#
