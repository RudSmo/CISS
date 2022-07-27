module Helper
using LinearAlgebra
export nkron, Idkron, fop, nfop, bop, HOBasis, ElNopUp, ElNopDw, PhNop

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
    Fermionic annihilation operator on site i on a lattice with L sites.
    """
    return [0 1; 0 0.]
end
function nfop(i,s,L,b)
    sz = [1 0; 0 -1.]
    id = Matrix(I,b,b)
    if i == 1 
        if s == 1
            return nkron([fop(),id,Idkron(2*(L-1),b)])
        else
            return nkron([sz,fop(),Idkron(2*(L-1),b)])
        end 
    else
        if s == 1
            return nkron([nkron([sz for j in 1:2*(i-1)]),fop(),id,Idkron(2*(L-i),b)])
        else
            return nkron([nkron([sz for j in 1:2*(i-1)]),sz,fop(),Idkron(2*(L-i),b)])
        end
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
function HOBasis(k::Int64,M::Int64)::Array
    """
    Generates harmonic oscillator basis state for M levels, occupied at k.
    Input:
        k ::Int64, occupied level
        M ::Int64, total number of levels
    """
    st = zeros(M)
    st[k+1] = 1.0
    return st
end
function ElNopUp(i::Int64,L::Int64)::Array
    """
    Constructs spin-resolved number operator of site i
    """
    nd = fop()*fop()'
    id = Matrix(I,2,2)
    return nkron([Idkron(i-1,2),Idkron(i-1,2),nd,id,Idkron(L-i,2),Idkron(L-i,2)])
end
function ElNopDw(i::Int64,L::Int64)::Array
    """
    Constructs spin-resolved number operator of site i
    """
    nd = fop()*fop()'
    id = Matrix(I,2,2)
    return nkron([Idkron(i-1,2),Idkron(i-1,2),id,nd,Idkron(L-i,2),Idkron(L-i,2)])
end
function PhNop(M::Int64,L::Int64,i::Int64)::Array
    """
    Bosonic number operator of boson with M levels, on site i with L sites.
    Input:
        M ::Int64, bosonic levels
        L ::Int64, lattice sites
        i ::Int64, site at which 
    """
    b = bop(M)
    bd = b'
    n = bd*b
    return nkron([Idkron(i-1,M),n,Idkron(L-i,M)])
end


end