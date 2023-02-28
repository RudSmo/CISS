include("../../../src/exact_TDSE/hamiltonian.jl")
using .Hamiltonian
using PyPlot
using LinearAlgebra
using SparseArrays
using Combinatorics 
using StatsBase

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


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
function Fransson(num_sites::Int64,vibStates::Int64,eps::Float64,t0::Float64,tPh::Float64,
                  wPh::Float64)::Array
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
    num = num_sites 
    dim = d^(num)*vibStates^(num-1)
    IdS = Matrix(I,d^(num),d^(num))
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
	cd_i = nFermi(cd,i,num)
	c_i = nFermi(c,i,num)

        hd = eps * OperatorProduct(cd_i,c_i)
        HS += kron(hd,IdB)
    end
    # Hopping
    for i in 1:num-1
	cd_i = nFermi(cd,i,num)
	c_i = nFermi(c,i+1,num)
        ht = -t0* OperatorProduct(cd_i,c_i)
	h_hop = ht+ht'
        HS += kron(h_hop,IdB)
	qv   = nkron(nBose(bd+b,i,num-1,vibStates))
	HSB -= tPh*kron(h_hop,qv)	
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
  	A = [("0,"^(sites-1)*"0",nkron([[1,0] for i in 1:sites]))]
	for i in 1:sites-1
      		on = ones(i)
      		ze = zeros(sites-i)
      		b = collect(Iterators.flatten([ze,on]))
	      	for bi in collect(myperms(b))[end:-1:1]
	  		str = ""
			for j in 1:sites
				k1 = bi[j]
				if k1 == 0
					str *= "0"
				elseif k1 == 1
					str *= "1"
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
  	push!(A,("1,"^(sites-1)*"1",nkron([[0,1] for i in 1:sites])))
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
	if iter == 1 && iter == N-1
		s0 = "|"
		s1 = "⟩|"
		s2 = "⟩"
	elseif iter == 1 && iter != N-1
		s0 = "|"
		s1 = "⟩|"
		s2 = ""
	elseif iter == N-1 && iter != 1 
		s0 = ""
		s1 = ""
                s2 = "⟩"
	end
	for a in A
		for b in B
			push!(conf,(s0*a[1]*s1*b[1]*s2,kron(a[2],b[2])))
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
	el_basis = ElectronBasis(nSites)	
	push!(state_confs,el_basis)
	for i in 1:nSites-1
		push!(state_confs,VibBasis(vibStates))
	end
	return MultiListKron(state_confs)
end


N = 2
vibStates = 160
T0 = 2.0 
t0 = T0
eps = 0.0
tPh = 0.2
wPh = 0.1


H = Fransson(N,vibStates,eps,t0,tPh,wPh)
#fil = open("Hamiltonian_Fransson_$(num_laps)x$(num_sites)_eps_$(E)_w0_$(wPh)_t0_$(T0)_t1_$(tPh)_tso_$(λ0)_tsoPh_$(λPh)_nv$(vibStates).txt","w")
#sH = sparse(H)
#println(sH)
#write(fil,H)


state_confs = GenerateFullBasis(N,vibStates)


ev = eigvecs(H)
println(size(ev))
GS = ev[:,1]

overlaps = []
confstr = []
for conf in state_confs
	cstr,cstate = conf
	ol = cstate' * GS
	if abs(ol) > 1e-3
		push!(overlaps,real.(ol))
		push!(confstr,cstr)
	end
end

so = sortperm(overlaps)
overlaps = overlaps[so][end:-1:1]
confstr = confstr[so][end:-1:1]

#for (cstr,ol) in zip(confstr,overlaps)
#	println(cstr," ",ol)
#end

phOccs = []
for cstr in confstr
	if string(cstr[9]) == "⟩"
		push!(phOccs,parse(Float64,string(cstr[8])))
	elseif length(string(cstr[8:end])) == 3
		push!(phOccs,parse(Float64,string(cstr[8:9])))
	else
		push!(phOccs,parse(Float64,string(cstr[8:10])))
	end
end




px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

ax.plot(phOccs,overlaps,linewidth=0,marker="o")
ax.set_ylabel("weight of configuration "*L"$\mathrm{Re}(a_n)$")
ax.set_xlabel("vibrational quantum number "*L"$n$")
plt.savefig("GS_nv_L2_Peierls_SimpleParams.pdf")




