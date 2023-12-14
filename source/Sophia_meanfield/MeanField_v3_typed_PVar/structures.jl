mutable struct HamiltonianParams


    onsight_e::Array{Float64}
    thopping::Array{Float64}
    lambda::Float64
    t_so::Float64
    sites_per_turn::Integer
    theta::Float64
    g_coupling::Float64 # array with two arrays ordered: j=1,j=2 then enumerating over n

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
    gPh         :: Array
    ws          :: Array
    X           :: Array # wie ist der displacement array sortiert??
    P           :: Array
    GL          :: Array
    GR          :: Array 
end

