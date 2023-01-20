module CustomStructures
export HamiltonianParams

mutable struct HamiltonianParams
    num_laps    :: Int64
    num_ions    :: Int64
    N           :: Int64
    Hamiltonian :: Function
    HamiltonStr :: String
    OnsiteEnergies :: Float64
    Hopping     :: Float64
    HoppingPh   :: Float64
    gPh         :: Float64
    ws          :: Array
    X           :: Array
    P           :: Array
    GL          :: Float64
    GR          :: Float64
end

end
