module Hamiltonian

include("helper.jl")

using .Helper

using LinearAlgebra
export Fransson,FranssonSymmetric,Holstein,HolsteinQuadratic,NNNTest,Hopping_t0t1,KondoImpurity,KondoImpurityNI


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
    dim = d^(2*num_sites)*vibStates^num_sites
    IdS = Matrix(I,d^(2*num_sites),d^(2*num_sites))
    IdB = Matrix(I,vibStates^num_sites,vibStates^num_sites)
    id = Matrix(I,2,2) 
    HS = zeros(dim,dim)*1im
    HB = zeros(dim,dim)   
    HSB = zeros(dim,dim)
    c = fop()'
    cd = c'
    nd = cd*c
    # On-site
    for i in 1:num_sites
        hu = eu * nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'
        hd = ed * nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)' 
        HS += kron(hu + hd,IdB)
    end
    # Hopping
    for i in 1:num_sites-1
        hu = t0* nfop(i,1,num_sites,d)*nfop(i+1,1,num_sites,d)'
        hd = t0* nfop(i,2,num_sites,d)*nfop(i+1,2,num_sites,d)'
        HS += kron(hu+hu'+hd+hd',IdB)
    end
    # SOC
    for i in 1:num_sites-2
        vs  = VSVec(num_sites,r_helix,c_helix,i,num_laps) 
        huu = 1im*λ0*vs[1,1]*nfop(i,1,num_sites,d)*nfop(i+2,1,num_sites,d)'
        hud = 1im*λ0*vs[1,2]*nfop(i,1,num_sites,d)*nfop(i+2,2,num_sites,d)'
        hdu = 1im*λ0*vs[2,1]*nfop(i,2,num_sites,d)*nfop(i+2,1,num_sites,d)'
        hdd = 1im*λ0*vs[2,2]*nfop(i,2,num_sites,d)*nfop(i+2,2,num_sites,d)'
        HS += kron(huu + hud + hdu + hdd + huu' + hud' + hdu' + hdd',IdB)
    end
    # Harmonic oscillators
    b = bop(vibStates)
    bd = b'
    bdb = bd*b
    for i in 1:num_sites
        HB += wPh*nkron([IdS,Idkron(i-1,vibStates),bdb,Idkron(num_sites-i,vibStates)]) 
    end
    # Electron-Phonon (NN)
    for i in 1:num_sites-1
        qv   = nkron([Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
	el_hop = nfop(i,1,num_sites,d)*nfop(i+1,1,num_sites,d)'
	el_hop += el_hop'
        HSB += tPh*kron(el_hop,qv)
	el_hop = nfop(i,2,num_sites,d)*nfop(i+1,2,num_sites,d)'
	el_hop += el_hop'
        HSB += tPh*kron(el_hop,qv)
    end
    # Electron-Phonon (SOC)
    for i in 1:num_sites-2
        vs  = VSVec(num_sites,r_helix,c_helix,i,num_laps) 
        qv  = nkron([Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
        huu = 1im*λPh*vs[1,1]*(nfop(i,1,num_sites,d)*nfop(i+2,1,num_sites,d))
        hud = 1im*λPh*vs[1,2]*(nfop(i,1,num_sites,d)*nfop(i+2,2,num_sites,d))
        hdu = 1im*λPh*vs[2,1]*(nfop(i,2,num_sites,d)*nfop(i+2,1,num_sites,d))
        hdd = 1im*λPh*vs[2,2]*(nfop(i,2,num_sites,d)*nfop(i+2,2,num_sites,d))
        HSB += kron(huu + hud + hdu + hdd + huu' + hud' + hdu' + hdd',qv) 
    end
    H = HS+HB+HSB
    return H
end
function FranssonSymmetric(num_sites::Int64,num_laps::Int64,vibStates::Int64,eps::Array,t0::Float64,tPh::Float64,
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
    dim = d^(2*num_sites)*vibStates^num_sites
    IdS = Matrix(I,d^(2*num_sites),d^(2*num_sites))
    IdB = Matrix(I,vibStates^num_sites,vibStates^num_sites)
    id = Matrix(I,2,2) 
    HS = zeros(dim,dim)*1im
    HB = zeros(dim,dim)   
    HSB = zeros(dim,dim)
    c = fop()'
    cd = c'
    nd = cd*c
    # On-site
    for i in 1:num_sites
        hu = eu * nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'
        hd = ed * nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)' 
        HS += kron(hu + hd,IdB)
    end
    # Hopping
    for i in 1:num_sites-1
        hu = -t0* nfop(i,1,num_sites,d)*nfop(i+1,1,num_sites,d)'
        hd = -t0* nfop(i,2,num_sites,d)*nfop(i+1,2,num_sites,d)'
        HS += kron(hu+hu'+hd+hd',IdB)
    end
    # SOC
    for i in 1:num_sites-2
        vs  = VSVec(num_sites,r_helix,c_helix,i,num_laps) 
        huu = -1im*λ0*vs[1,1]*nfop(i,1,num_sites,d)*nfop(i+2,1,num_sites,d)'
        hud = -1im*λ0*vs[1,2]*nfop(i,1,num_sites,d)*nfop(i+2,2,num_sites,d)'
        hdu = -1im*λ0*vs[2,1]*nfop(i,2,num_sites,d)*nfop(i+2,1,num_sites,d)'
        hdd = -1im*λ0*vs[2,2]*nfop(i,2,num_sites,d)*nfop(i+2,2,num_sites,d)'
        HS += kron(huu + hud + hdu + hdd + huu' + hud' + hdu' + hdd',IdB)
    end
    # Harmonic oscillators
    b = bop(vibStates)
    bd = b'
    bdb = bd*b
    for i in 1:num_sites
        HB += wPh*nkron([IdS,Idkron(i-1,vibStates),bdb,Idkron(num_sites-i,vibStates)]) 
    end
    # Electron-Phonon (NN)
    for i in 1:num_sites-1
        qv   = nkron([Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
	el_hop = nfop(i,1,num_sites,d)*nfop(i+1,1,num_sites,d)'
	el_hop += el_hop'
        HSB -= tPh*kron(el_hop,qv)
	el_hop = nfop(i,2,num_sites,d)*nfop(i+1,2,num_sites,d)'
	el_hop += el_hop'
        HSB -= tPh*kron(el_hop,qv)
    end
    qv   = nkron([Idkron(num_sites-1,vibStates),bd+b])
    el_hop = nfop(num_sites,1,num_sites,d) * nfop(num_sites-1,1,num_sites,d)'
    el_hop += el_hop'
    HSB = -tPh*kron(el_hop,qv)
    el_hop = nfop(num_sites,2,num_sites,d) * nfop(num_sites-1,2,num_sites,d)'
    el_hop += el_hop'
    HSB = -tPh*kron(el_hop,qv)


    # Electron-Phonon (SOC)
    for i in 1:num_sites-2
        vs  = VSVec(num_sites,r_helix,c_helix,i,num_laps) 
        qv  = nkron([Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
        huu = -1im*λPh*vs[1,1]*(nfop(i,1,num_sites,d)*nfop(i+2,1,num_sites,d))
        hud = -1im*λPh*vs[1,2]*(nfop(i,1,num_sites,d)*nfop(i+2,2,num_sites,d))
        hdu = -1im*λPh*vs[2,1]*(nfop(i,2,num_sites,d)*nfop(i+2,1,num_sites,d))
        hdd = -1im*λPh*vs[2,2]*(nfop(i,2,num_sites,d)*nfop(i+2,2,num_sites,d))
        HSB += kron(huu + hud + hdu + hdd + huu' + hud' + hdu' + hdd',qv) 
    end
    H = HS+HB+HSB
    return H
end


function Holstein(num_sites::Int64,vibStates::Int64,eps::Array,t0::Float64,gPh::Float64,wPh::Float64)::Array
    """
    Input:
        num_sites:: Int64, lattice sites
        vibStates:: Int64, bosonic levels
        eps      :: Array, on-site energies
        t0       :: Float64, nearest-neighbor hopping
        gph      :: Electron-phonon coupling
        wPh      :: Harmonic oscillator frequency
    """
    d = 2
    eu,ed = eps
    dim = d^(2*num_sites)*vibStates^num_sites
    IdS = Matrix(I,d^(2*num_sites),d^(2*num_sites))
    IdB = Matrix(I,vibStates^num_sites,vibStates^num_sites)
    id = Matrix(I,2,2) 
    HS = zeros(dim,dim)*1im
    HB = zeros(dim,dim)   
    HSB = zeros(dim,dim)
    c = fop()'
    cd = c'
    nd = cd*c
    # On-site
    for i in 1:num_sites
        hu = eu * nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'
        hd = ed * nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)' 
        HS += kron(hu + hd,IdB)
    end
    # Hopping
    for i in 1:num_sites-1
        hu = -t0* nfop(i,1,num_sites,d)*nfop(i+1,1,num_sites,d)'
        hd = -t0* nfop(i,2,num_sites,d)*nfop(i+1,2,num_sites,d)'
        HS += kron(hu+hu'+hd+hd',IdB)
    end
    # Harmonic oscillators
    b = bop(vibStates)
    bd = b'
    bdb = bd*b
    for i in 1:num_sites
        HB += wPh*nkron([IdS,Idkron(i-1,vibStates),bdb,Idkron(num_sites-i,vibStates)]) 
    end
    # Electron-Phonon (NN)
    for i in 1:num_sites-1
        qv   = nkron([Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
        HSB += gPh*kron(nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'+nfop(i+1,1,num_sites,d)*nfop(i,1,num_sites,d)',qv)
        HSB += gPh*kron(nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)'+nfop(i+1,2,num_sites,d)*nfop(i,2,num_sites,d)',qv)
    end
    return HS+HB+HSB
end

function HolsteinQuadratic(num_sites::Int64,vibStates::Int64,eps::Array,t0::Float64,gPh::Float64,gPh2::Float64,wPh::Float64)::Array
    """
    Input:
        num_sites:: Int64, lattice sites
        vibStates:: Int64, bosonic levels
        eps      :: Array, on-site energies
        t0       :: Float64, nearest-neighbor hopping
        gph      :: Linear electron-phonon coupling
        gph2     :: Quadratic electron-phonon coupling
        wPh      :: Harmonic oscillator frequency
    """
    d = 2
    eu,ed = eps
    dim = d^(2*num_sites)*vibStates^num_sites
    IdS = Matrix(I,d^(2*num_sites),d^(2*num_sites))
    IdB = Matrix(I,vibStates^num_sites,vibStates^num_sites)
    id = Matrix(I,2,2) 
    HS = zeros(dim,dim)*1im
    HB = zeros(dim,dim)   
    HSB = zeros(dim,dim)
    c = fop()'
    cd = c'
    nd = cd*c
    # On-site
    for i in 1:num_sites
        hu = eu * nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'
        hd = ed * nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)' 
        HS += kron(hu + hd,IdB)
    end
    # Hopping
    for i in 1:num_sites-1
        hu = -t0* nfop(i,1,num_sites,d)*nfop(i+1,1,num_sites,d)'
        hd = -t0* nfop(i,2,num_sites,d)*nfop(i+1,2,num_sites,d)'
        HS += kron(hu+hu'+hd+hd',IdB)
    end
    # Harmonic oscillators
    b = bop(vibStates)
    bd = b'
    bdb = bd*b
    for i in 1:num_sites
        HB += wPh*nkron([IdS,Idkron(i-1,vibStates),bdb,Idkron(num_sites-i,vibStates)]) 
    end
    # Electron-Phonon (NN)
    for i in 1:num_sites-1
        qv   = nkron([Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
        HSB += gPh*kron(nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'+nfop(i+1,1,num_sites,d)*nfop(i,1,num_sites,d)',qv)
        HSB += gPh*kron(nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)'+nfop(i+1,2,num_sites,d)*nfop(i,2,num_sites,d)',qv)
    end
    # Electron-Phonon (NN)
    for i in 1:num_sites-1
        qv   = nkron([Idkron(i-1,vibStates),(bd+b)*(bd+b),Idkron(num_sites-i,vibStates)])
        HSB += gPh2*kron(nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'+nfop(i+1,1,num_sites,d)*nfop(i,1,num_sites,d)',qv)
        HSB += gPh2*kron(nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)'+nfop(i+1,2,num_sites,d)*nfop(i,2,num_sites,d)',qv)
    end 
    return HS+HB+HSB
end

function NNNTest(num_sites,eps,t1,t2,vibStates)
    eu,ed = eps
    t2uu,t2ud,t2dd = t2
    d = 2
    dim = d^(2*num_sites)*vibStates^num_sites
    IdS = Matrix(I,d^(2*num_sites),d^(2*num_sites))
    IdB = Matrix(I,vibStates^num_sites,vibStates^num_sites)
    id = Matrix(I,2,2) 
    HS = zeros(dim,dim)*1im
    c = fop()'
    cd = c'
    nd = cd*c
    # On-site
    for i in 1:num_sites
        hu = eu * nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'
        hd = ed * nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)' 
        HS += kron(hu + hd,IdB)
    end
    # Hopping
    for i in 1:num_sites-1
        hu = -t1* nfop(i,1,num_sites,d)*nfop(i+1,1,num_sites,d)'
        hd = -t1* nfop(i,2,num_sites,d)*nfop(i+1,2,num_sites,d)'
        HS += kron(hu+hu'+hd+hd',IdB)
    end
    # NNN hopping
    for i in 1:num_sites-2
        huu = -t2uu* nfop(i,1,num_sites,d)*nfop(i+2,1,num_sites,d)'
        hud = -t2ud* nfop(i,1,num_sites,d)*nfop(i+2,2,num_sites,d)'
        hdu = -conj(t2ud)* nfop(i,2,num_sites,d)*nfop(i+2,1,num_sites,d)'
        hdd = -t2dd* nfop(i,2,num_sites,d)*nfop(i+2,2,num_sites,d)'
         

        HS += kron(huu+huu'+hud+hud'+hdu+hdu'+hdd+hdd',IdB)
    end
    
    return HS
end
function Hopping_t0t1(num_sites,eps,t0,t1,vibStates)
    eu,ed = eps
    d = 2
    dim = d^(2*num_sites)*vibStates^num_sites
    IdS = Matrix(I,d^(2*num_sites),d^(2*num_sites))
    IdB = Matrix(I,vibStates^num_sites,vibStates^num_sites)
    id = Matrix(I,2,2)
    HS = zeros(dim,dim)*1im
    c = fop()'
    cd = c'
    nd = cd*c
    # On-site
    for i in 1:num_sites
        hu = eu * nfop(i,1,num_sites,d)*nfop(i,1,num_sites,d)'
        hd = ed * nfop(i,2,num_sites,d)*nfop(i,2,num_sites,d)'
        HS += kron(hu + hd,IdB)
    end
    # Hopping
    for i in 1:num_sites-1
        hu = -t0* nfop(i,1,num_sites,d)*nfop(i+1,1,num_sites,d)'
        hd = -t0* nfop(i,2,num_sites,d)*nfop(i+1,2,num_sites,d)'

	hud = -t1 * nfop(i,1,num_sites,d)*nfop(i+1,2,num_sites,d)'
	hdu = -t1 * nfop(i,2,num_sites,d)*nfop(i+1,1,num_sites,d)'
        HS += kron(hu+hu'+hd+hd'+hud+hud'+hdu+hdu',IdB)
    end

    return HS
end

function KondoImpurity(L::Int64,t0::Float64,U::Float64,JK::Float64)
    dim = 4^L
    d = 2
    HS = zeros(dim,dim) * (1+0.0*1im)
    IS = Matrix(I,dim,dim)
    IB = Matrix(I,4,4)
    # Hopping in the chain
    for i in 1:L-1
        He = -t0 * nfop(i,1,L,d) * nfop(i+1,1,L,d)' 
        He += -t0 * nfop(i,2,L,d) * nfop(i+1,2,L,d)' 
        HS += He + He' 
    end
    HT = kron(IB,HS)
    # Coulomb on-site in the chain
    for i in 1:L
        HS += U*nfop(i,1,L,d)*nfop(i,1,L,d)' * nfop(i,2,L,d)*nfop(i,2,L,d)' 
    end
    HT += kron(IB,HS)


    Ddag_up = nfop(1,1,1,d)
    Ddag_dn = nfop(1,2,1,d)
    D_up = Ddag_up'
    D_dn = Ddag_dn' 
    Cdag_up = nfop(1,1,L,d)
    Cdag_dn = nfop(1,2,L,d)
    C_up = Cdag_up'
    C_dn = Cdag_dn'
     
    HT += JK * (1/2 * ( kron( Ddag_dn*D_up , Cdag_up*C_dn ) + kron( Ddag_up*D_dn , Cdag_dn*C_up ) )
               +1/4 * ( kron( Ddag_up*D_up - Ddag_dn*D_dn , Cdag_up*C_up - Cdag_dn*C_dn ) ) )


end

function KondoImpurityNI(L::Int64,t0::Float64,JK::Float64)
    dim = 2*L
    d = 2
    HS = zeros(dim,dim) * (1+0.0*1im)
    IS = Matrix(I,dim,dim)
    IB = Matrix(I,4,4)
    # Hopping in the chain
    for i in 1:2*(L-1)
        HS[i,i+2] = -t0
        HS[i+2,i] = -t0
    end
    HT = kron(IB,HS)
    

    Ddag_up = nfop(1,1,1,d)
    Ddag_dn = nfop(1,2,1,d)
    D_up = Ddag_up'
    D_dn = Ddag_dn'
    
    Htmp = zeros(dim,dim)
    Hud = Htmp
    Hdu = Htmp
    Huu = Htmp
    Hdd = Htmp
    
    Hud[1,2] = 1.0
    Hdu[2,1] = 1.0
    Huu[1,1] = 1.0
    Hdd[2,2] = 1.0
     
    HT += JK * (1/2 * ( kron( Ddag_dn*D_up , Hdu ) + kron( Ddag_up*D_dn , Hud ) )
               +1/4 * ( kron( Ddag_up*D_up - Ddag_dn*D_dn , Huu - Hdd ) ) )


end



end
