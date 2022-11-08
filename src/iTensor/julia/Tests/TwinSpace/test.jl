using ITensors
using LinearAlgebra
using PyPlot
using Einsum


function nkron(A::Array)::Any
    """
    Tensor product of all arrays, stored in A.
    Input:
        A ::Array = [A1,A2,...,An]
    Output:
        C ::Array = A1⊗A2⊗...⊗An
    """
    if length(A) == 1
        return A[1]
    else
        return kron(A[1],nkron(A[2:end]))
    end
end

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
function fop()::Array
    """
    Single fermionic annihilation operator.
    """
    return [0 1; 0 0.]
end
function nfop_twinSpace(i,s,L)
    """
    Fermionic annihilation operator on site i and spin s on a lattice with L sites and dimension b.
    """
    sz = [1 0; 0 -1.]
    id = Matrix(I,2,2)
    if i == 1 
        if s == 1
            return nkron([fop(),id,Idkron(2*(L-1),2)])
        else
            return nkron([id,fop(),Idkron(2*(L-1),2)])
        end 
    else #if L > i > 1
        if s == 1
            id_list_L = []
            for k in 1:i-1
                push!(id_list_L,sz)
                push!(id_list_L,id)
            end
            return nkron([nkron(id_list_L),fop(),id,Idkron(2*(L-i),2)])
        else
            id_list_R = []
            for k in 1:i-1
                push!(id_list_R,id)
                push!(id_list_R,sz)
            end
            return nkron([nkron(id_list_R),id,fop(),Idkron(2*(L-i),2)])
        end
    end
end


L = 3
i = 1
j = 3

H_direct = nfop_twinSpace(i,1,L)' * nfop_twinSpace(j,1,L) - nfop_twinSpace(i,2,L)' * nfop_twinSpace(j,2,L) 
H_direct = H_direct + H_direct'
fig, ax = plt.subplots(2,1)
ax[1].imshow(real.(H_direct))

function ITensor_Operator_twinSpace(op_p::String,i::Int64,s::Int64,L::Int64)::Prod{Op}
    A = Prod{Op}()
    k = 1
    for j in 1:i-1
        if s == 1 # physical space
            A *= "F",k
            A *= "Id",k+1
            k += 2
        else
            A *= "Id",k
            A *= "F",k+1
            k += 2 
        end
    end
    if s == 1
        A *= "$(op_p)",k
        A *= "Id",k+1
        k += 2
    else
        A *= "Id",k
        A *= "$(op_p)",k+1
        k += 2 
    end
    for j in i+1:L
        A *= "Id",k
        A *= "Id",k+1
        k += 2
    end
    return A  
end
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
    return A[end:-1:1]  
end


function ITensor_Prod(ops,k)
    if length(ops) == 1
        return ITensors.terms(ops[1])[k]
    else
        return ITensors.terms(ops[1])[k] * ITensor_Prod(ops[2:end],k)
    end
end

function ITensor_ProdString(ops,k)
    new_op = ops[1][k]
    n = length(ops)
    for j in 2:n
        new_op *= " * $(ops[j][k])"
    end
    return new_op
end

function ITensor_Prod_twinSpace(ops::Vector{Prod{Op}},sites::Vector{Index{Int64}})::MPO
    m = length(ops[1])   
    M = Vector{ITensor}()
    for j in 1:m
        Mp = MPO(ITensor_Prod(ops,j),sites)[j]
        push!(M,Mp)
    end
    return MPO(M)        
end

function ITensor_ProdString_twinSpace(ops::Vector{Vector{String}},sites::Vector{Index{Int64}})::MPO
    
    m = length(ops[1])
#    str_op = [ITensor_ProdString(ops,j) for j in 1:m]
#    return ITensors.MPO(sites,str_op)
    op_fin = Vector{ITensor}()
    for j in 1:m
        str_tmp = ITensor_ProdString(ops,j)
        op_tmp = ITensors.op(str_tmp,sites[j])
        push!(op_fin,op_tmp)
    end
    return ITensors.MPO(op_fin)
end


ITensors.op(::OpName"Id",::SiteType"Fermion") = [1 0; 0 1.]
sites = ITensors.siteinds("Fermion",2*L)
cdag_ts_p = ITensor_OperatorString_twinSpace("Cdag",i,1,L)
c_ts_p = ITensor_OperatorString_twinSpace("C",j,1,L)
ops_p = [cdag_ts_p,c_ts_p]
mpo_phys = ITensor_ProdString_twinSpace(ops_p,sites)

cdag_ts_p = ITensor_OperatorString_twinSpace("Cdag",j,1,L)
c_ts_p = ITensor_OperatorString_twinSpace("C",i,1,L)
ops_p = [cdag_ts_p,c_ts_p]
mpo_phys = mpo_phys + ITensor_ProdString_twinSpace(ops_p,sites)


cdag_ts_a = ITensor_OperatorString_twinSpace("Cdag",i,2,L)
c_ts_a = ITensor_OperatorString_twinSpace("C",j,2,L)
ops_a = [cdag_ts_a,c_ts_a]
mpo_ancilla = ITensor_ProdString_twinSpace(ops_a,sites)

cdag_ts_a = ITensor_OperatorString_twinSpace("Cdag",j,2,L)
c_ts_a = ITensor_OperatorString_twinSpace("C",i,2,L)
ops_a = [cdag_ts_a,c_ts_a]
mpo_ancilla = mpo_ancilla + ITensor_ProdString_twinSpace(ops_a,sites) 

H_mpo = mpo_phys - mpo_ancilla

Hmpo1 = array(H_mpo[1]); 
Hmpo2 = array(H_mpo[2]); 
Hmpo3 = array(H_mpo[3]); 
Hmpo4 = array(H_mpo[4]); 
Hmpo5 = array(H_mpo[5]); 
Hmpo6 = array(H_mpo[6])
@einsum HT[i,j,l,m,o,p,r,s,u,v,x,y] := Hmpo1[i,j,k] * Hmpo2[l,m,n,k] * Hmpo3[o,p,q,n] * Hmpo4[r,s,t,q] * Hmpo5[u,v,w,t] * Hmpo6[x,y,w]

#HT = permutedims(HT,(11,9,7,5,3,1,12,10,8,6,4,2))
HT = permutedims(HT,(1,3,5,7,9,11,2,4,6,8,10,12))
#!!!yields correct hamiltonian!!!
HT = reshape(HT,size(H_direct))

ax[2].imshow(real.(HT))
plt.show()

