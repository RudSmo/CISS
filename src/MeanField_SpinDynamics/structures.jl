module CustomStructures
export HamiltonianParams

mutable struct HamiltonianParams
    N           :: Int64
    Hamiltonian :: Function
    HamiltonStr :: String
    OnsiteEnergies :: Array
    Hopping     :: Float64
    Jsd 	:: Float64
    mag         :: Array
    pos         :: Array
    Hext        :: Array
    tL          :: Array
    tR          :: Array 
end

end
