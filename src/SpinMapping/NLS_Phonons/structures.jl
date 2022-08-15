module Structures

export HamiltonianParameters, PhaseSpaceVariables, SimulationParameters

struct HamiltonianParameters
    """
    Hamiltonian = Electronic Hamiltonian
    gamma       = Electron-phonon couplings
    omega       = Phonon energies
    Electron-phonon coupling to be assumed as local Holstein
    """
    Hamiltonian ::Array{Float64}
    gamma       ::Array{Float64}
    omega       ::Array{Float64}
    beta        ::Float64
end

mutable struct PhaseSpaceVariables
    """
    Ordering (Example): Three-level system
    X = Mapping coordinate (electrons)
    P = Mapping momentum (electrons)
    x = Phonon coordinates 
    p = Phonon momenta
    """
    X       ::Array{ComplexF64}
    P       ::Array{ComplexF64}
    x       ::Array{ComplexF64}
    p       ::Array{ComplexF64}
    gs      ::Float64
end

struct SimulationParameters
    Nlvl    ::Int64
    sphere  ::String
    ntraj   ::Int64
    dt      ::Float64
    Nstp    ::Int64
    Nsave   ::Int64 
    Nflush  ::Int64
    fnE     ::String
    fnP     ::String
end

end