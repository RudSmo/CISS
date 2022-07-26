module Hamiltonian
#this is a tiny change

include("helper.jl")

using .Helper

using LinearAlgebra
export Fransson


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
                  λ0::Float64,λPh::Float64,r_helix::Float64,c_helix::Float64)::Array
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
        hu = eu * nkron([Idkron(i-1,d),Idkron(i-1,d),cd * c,id,Idkron(num_sites-i,d),Idkron(num_sites-i,d)])
        hd = ed * nkron([Idkron(i-1,d),Idkron(i-1,d),id,cd * c,Idkron(num_sites-i,d),Idkron(num_sites-i,d)])
        HS += kron(hu + hd,IdB)
    end
    # Hopping
    for i in 1:num_sites-1
        hu = -t0*nkron([Idkron(i-1,d),Idkron(i-1,d),cd,id,c,id,Idkron(num_sites-1-i,d),Idkron(num_sites-1-i,d)]) 
        hd = -t0*nkron([Idkron(i-1,d),Idkron(i-1,d),id,cd,id,c,Idkron(num_sites-1-i,d),Idkron(num_sites-1-i,d)]) 
        HS += kron(hu+hu'+hd+hd',IdB)
    end
    # SOC
    for i in 1:num_sites-2
        vs  = VSVec(num_sites,r_helix,c_helix,i,num_laps) 
        huu = -1im*λ0*vs[1,1]*nkron([Idkron(i-1,d),Idkron(i-1,d),cd,id, id,id, c,id, Idkron(num_sites-2-i,d),Idkron(num_sites-2-i,d)]) 
        hud = -1im*λ0*vs[1,2]*nkron([Idkron(i-1,d),Idkron(i-1,d),cd,id, id,id, id,c, Idkron(num_sites-2-i,d),Idkron(num_sites-2-i,d)]) 
        hdu = -1im*λ0*vs[2,1]*nkron([Idkron(i-1,d),Idkron(i-1,d),id,cd, id,id, c,id, Idkron(num_sites-2-i,d),Idkron(num_sites-2-i,d)]) 
        hdd = -1im*λ0*vs[2,2]*nkron([Idkron(i-1,d),Idkron(i-1,d),id,cd, id,id, id,c, Idkron(num_sites-2-i,d),Idkron(num_sites-2-i,d)]) 
        HS += kron(huu + hud + hdu + hdd + huu' + hud' + hdu' + hdd',IdB)
    end
    # Harmonic oscillators
    b = bop(vibStates)
    bd = b'
    bdb = bd*b
    for i in 1:num_sites
        HB += nkron([IdS,Idkron(i-1,vibStates),bdb,Idkron(num_sites-i,vibStates)]) 
    end
    # Electron-Phonon (NN)
    for i in 1:num_sites-1
        HSB -= tPh*nkron([Idkron(i-1,d),Idkron(i-1,d),cd,id,c,id,Idkron(num_sites-i-1,d),Idkron(num_sites-i-1,d),Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
        HSB -= tPh*nkron([Idkron(i-1,d),Idkron(i-1,d),id,cd,id,c,Idkron(num_sites-i-1,d),Idkron(num_sites-i-1,d),Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])     
    end
    # Electron-Phonon (SOC)
    for i in 1:num_sites-2
        vs   = VSVec(num_sites,r_helix,c_helix,i,num_laps) 

        huu = -1im*λPh*vs[1,1]*nkron([Idkron(i-1,d),Idkron(i-1,d),cd,id, id,id, c,id,Idkron(num_sites-i-2,d),Idkron(num_sites-i-2,d),Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
        hud = -1im*λPh*vs[1,2]*nkron([Idkron(i-1,d),Idkron(i-1,d),cd,id, id,id, id,c,Idkron(num_sites-i-2,d),Idkron(num_sites-i-2,d),Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])     
        hdu = -1im*λPh*vs[2,1]*nkron([Idkron(i-1,d),Idkron(i-1,d),id,cd, id,id, c,id,Idkron(num_sites-i-2,d),Idkron(num_sites-i-2,d),Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])
        hdd = -1im*λPh*vs[2,2]*nkron([Idkron(i-1,d),Idkron(i-1,d),id,cd, id,id, id,c,Idkron(num_sites-i-2,d),Idkron(num_sites-i-2,d),Idkron(i-1,vibStates),bd+b,Idkron(num_sites-i,vibStates)])     
        HSB += huu + hud + hdu + hdd + huu' + hud' + hdu' + hdd' 
    end
    return HS+HB+HSB
end
end
