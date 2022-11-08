using ITensors
using LinearAlgebra
using ITensorTDVP
using KrylovKit
using LinearAlgebra
using Observers

ITensors.op(::OpName"Id",::SiteType"Fermion") = [1 0; 0 1.]

function ITensor_OperatorString_twinSpace(op_p::String,i::Int64,s::Int64,L::Int64)::Vector{String}
    A = Vector{String}()
    k = 1
    for j in 1:i-1
        if s == 1 # physical space
            push!(A,"F")
            push!(A,"Id")
        else
            push!(A,"Id")
            push!(A,"F")
        end
    end
    if s == 1
        push!(A,"$(op_p)")
        push!(A,"Id")
    else
        push!(A,"Id")
        push!(A,"$(op_p)")
    end
    for j in i+1:L
        push!(A,"Id")
        push!(A,"Id")
    end
    return A#[end:-1:1]  
end
function ITensor_ProdString(ops,k)
    new_op = ops[1][k]
    n = length(ops)
    for j in 2:n
        new_op *= " * $(ops[j][k])"
    end
    return new_op
end
function ITensor_ProdString_twinSpace(ops::Vector{Vector{String}},sites::Vector{Index{Int64}})::MPO
    
    m = length(ops[1])
    op_fin = Vector{ITensor}()
    for j in 1:m   
        str_tmp = ITensor_ProdString(ops,j)
        op_tmp = ITensors.op(str_tmp,sites[j]) 
        push!(op_fin,op_tmp)
    end
    return ITensors.MPO(op_fin)
end


N = 3
t0 = -1.0
dt = 0.005
time_stop = 10.0
sites = siteinds("Fermion",2*N)
state = []
for i in 1:N
	if i == 1
		push!(state,"Occ")
        push!(state,"Occ")
	else
		push!(state,"Emp")
        push!(state,"Emp")      
	end
end

function HamiltonianMPO_twinSpace(L,sites,t0)
    mpo_p = 0.0 * ITensors.MPO(sites,["Id" for i in 1:2*L]) 
    mpo_a = 0.0 * ITensors.MPO(sites,["Id" for i in 1:2*L])
    for i in 1:L-1
        
        cdag_p = ITensor_OperatorString_twinSpace("Cdag",i,1,L)
        cdag_a = ITensor_OperatorString_twinSpace("Cdag",i,2,L)
        c_p = ITensor_OperatorString_twinSpace("C",i+1,1,L)
        c_a = ITensor_OperatorString_twinSpace("C",i+1,2,L)

        ops_p = [cdag_p,c_p]
        mpo_p += t0 * ITensor_ProdString_twinSpace(ops_p,sites) 
        ops_a = [cdag_a,c_a]
        mpo_a += t0 * ITensor_ProdString_twinSpace(ops_a,sites) 

        cdag_p = ITensor_OperatorString_twinSpace("Cdag",i+1,1,L)
        cdag_a = ITensor_OperatorString_twinSpace("Cdag",i+1,2,L)

        c_p = ITensor_OperatorString_twinSpace("C",i,1,L)
        c_a = ITensor_OperatorString_twinSpace("C",i,2,L)

        ops_p = [cdag_p,c_p]
        mpo_p += t0 * ITensor_ProdString_twinSpace(ops_p,sites) 
        ops_a = [cdag_a,c_a]
        mpo_a += t0 * ITensor_ProdString_twinSpace(ops_a,sites)
    end
    return mpo_p - mpo_a
end


psi = productMPS(sites,state)
H   = HamiltonianMPO_twinSpace(N,sites,t0)

outputlevel = 1
nsite = 1
start_linkdim = 12
maxdim = 200
cutoff = 1e-9
tol = 1e-15
krylov_kwargs = (; tol=tol, eager=true)

function step(; sweep, bond, half_sweep)
    if bond == 1 && half_sweep == 2
      return sweep
    end
    return nothing
  end
function current_time(; current_time, bond, half_sweep)
    if bond == 1 && half_sweep == 2
      return current_time
    end
    return nothing
  end
function measure_nd(; psi, bond, half_sweep)
    if bond == 1 && half_sweep == 2
        Nd = []
        for j in 1:2:2*N
            push!(Nd,expect(psi, "N"; sites=j))
        end
        return Nd
    end
    return nothing
end

function return_state(; psi, bond, half_sweep)
    if bond == 1 && half_sweep == 2
      return psi
    end
    return nothing
  end
obs = Observer("steps" => step, "times" => current_time, "psis" => return_state, "Nd" => measure_nd)
psi_t = tdvp(
    H,
    time_stop*1im,
    psi;
    time_step=1im*dt,
    cutoff,
    normalize=false,
    maxdim,
    (observer!)=obs,
    outputlevel=1)    
res = results(obs)
steps = res["steps"]
times = res["times"]
psis = res["psis"]
Nd = res["Nd"]

fil = open("TDVP_L$(N)_dt$(dt).dat","w")
for i in 1:length(Nd)
    write(fil,"$(imag(times[i])) ")
    for j in 1:length(Nd[i])
        write(fil,"$(Nd[i][j]) ")
    end
    write(fil,"\n")
end


