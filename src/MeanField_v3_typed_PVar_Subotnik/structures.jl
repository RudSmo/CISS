mutable struct HamiltonianParams
    N           :: Int64
    Hamiltonian :: Function
    HamiltonStr :: String
    A		:: Float64
    B           :: Float64
    Delta	:: Float64
    kappa_x	:: Array
    chi_x	:: Array
    X		:: Array
    P           :: Array
    GL          :: Array
    GR          :: Array 
end

