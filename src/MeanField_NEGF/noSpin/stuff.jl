function G0_R(HS,E,eta)
    """
    Central system retarded Green function.
    """
    id = Matrix(I,size(H))
    return inv(E*id-HS-1im*eta*id)
end

function G0_A(G0_R)
    """
    Central system advanced Green function.
    """
    return G0_R'
end

function Sigma_R(Gamma)
    """
    Wide band limit.
    """
    return -1im/2 * Gamma
end

function Sigma_lsr(Gamma,E,T)
    """
    Gamma is a function of energy.
    """
    return 1im * Gamma * Fermi(E,T)
end

function Sigma_gtr(Gamma,E,T)
    return -1im * Gamma * (1-Fermi(E,T))
end

function Gamma_SI(Gamma,E,t,ec)
    """
    Level-width function for a semi-infinite chain.
    """
    if 4*abs(t)^2 > (E-ec)^2
        return Gamma^2/(abs(t)^2) * sqrt(4*abs(t)^2-(E-ec)^2)
    else
        return 0.0
    end
end
