module CustomStructures
export HamiltonianParams

mutable struct HamiltonianParams
    num_laps    :: Int64
    num_ions    :: Int64
    N           :: Int64
    Hamiltonian :: Function
    HamiltonStr :: String
    OnsiteEnergies :: Array
    a_hel       :: Float64
    c_hel       :: Float64
    Hopping     :: Float64
    HoppingPh   :: Array
    SOC         :: Float64
    SOCPh       :: Float64
    Kph         :: Float64
    Mph		:: Float64
    ws		:: Float64
    X           :: Array
    P           :: Array
    GL          :: Array
    GR          :: Array 
end

end
