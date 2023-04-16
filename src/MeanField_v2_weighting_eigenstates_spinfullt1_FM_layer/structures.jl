module CustomStructures
export HamiltonianParams

mutable struct HamiltonianParams
    fm_sites	:: Array
    mol_sites	:: Array
    num_laps    :: Int64
    num_ions    :: Int64
    N           :: Int64
    Hamiltonian :: Function
    HamiltonStr :: String
    OnsiteEnergiesFM :: Array
    OnsiteEnergies :: Array
    coup_site	:: Int64
    coup_en	:: Array
    a_hel       :: Float64
    c_hel       :: Float64
    Hopping_FM	:: Float64
    Hopping     :: Float64
    HoppingPh   :: Array
    SOC         :: Float64
    SOCPh       :: Float64
    gPh         :: Array
    ws          :: Array
    X           :: Array
    P           :: Array
    GL          :: Array
    GR          :: Array 
    mag		:: Array
    Jsd		:: Float64
end

end
