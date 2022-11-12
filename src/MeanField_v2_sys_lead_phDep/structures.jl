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
    HoppingPh   :: Float64
    SOC         :: Float64
    SOCPh       :: Float64
    gPh         :: Array
    ws          :: Array
    X           :: Array
    P           :: Array 
    coupL       :: Array
    coupR       :: Array
end

end
